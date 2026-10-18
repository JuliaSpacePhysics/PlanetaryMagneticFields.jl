coeffs := "data/coeffs"
lif := "mattkjames7/libinternalfield"
pmf := "AnkitBarik/planetMagFields"
# planetMagFields models absent from libinternalfield, as <planet>/<model>.
pmf_models := "mercury/wardinski2019 saturn/cassini11+"
# planetMagFields files duplicating libinternalfield models (some in other units or normalization).
pmf_dupes := "earth_igrf14 ganymede_kivelson2002 jupiter_jrm09 jupiter_jrm33 jupiter_vip4 mercury_anderson2012 mercury_thebault2018 neptune_connerney1991 neptune_holme1996 saturn_cassini11 saturn_cassinisoi uranus_connerney1987 uranus_herbert2009 uranus_holme1996"
# Local files ahead of libinternalfield (IGRF-14); kept as-is by check and sync.
pinned := "earth/igrf2020.dat earth/igrf2025.dat"

# Report differences between upstream coefficients and data/coeffs
check: (_sync "check")

# Update data/coeffs from upstream
sync: (_sync "apply")

_sync mode:
    #!/usr/bin/env bash
    set -euo pipefail
    up=$(mktemp -d); trap 'rm -rf "$up"' EXIT

    git clone -q --depth 1 --filter=blob:none --sparse https://github.com/{{lif}} "$up/.lif"
    git -C "$up/.lif" sparse-checkout set data/coeffs
    cp -R "$up/.lif/data/coeffs/." "$up/"
    rm -rf "$up/.lif"

    raw="https://raw.githubusercontent.com/{{pmf}}/HEAD/planetmagfields/data"
    for m in {{pmf_models}}; do
        curl -sfL "$raw/${m%%/*}_${m#*/}.dat" | awk 'NF{print $1"\t"$2"\t"$3"\t"$4}' > "$up/$m.dat"
    done
    known=" {{pmf_dupes}} $(echo {{pmf_models}} | tr / _) "
    gh api "repos/{{pmf}}/contents/planetmagfields/data" --jq '.[].name' | while read -r f; do
        [[ $known == *" ${f%.dat} "* ]] || echo "unreviewed planetMagFields model: $f" >&2
    done

    for f in {{pinned}}; do cp "{{coeffs}}/$f" "$up/$f"; done

    if [[ {{mode}} == check ]]; then
        diff -rq "$up" {{coeffs}} && echo "{{coeffs}} is up to date"
    else
        rsync -rc --out-format='updated %n' "$up/" {{coeffs}}/
    fi
