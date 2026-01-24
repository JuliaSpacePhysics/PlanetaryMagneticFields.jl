planetMagFields_data:
    #!/bin/bash
    set -e
    # Create data directory if it doesn't exist
    mkdir -p data
    # Download all data files
    base_url="https://raw.githubusercontent.com/AnkitBarik/planetMagFields/master/planetmagfields/data"
    # Fallback: download known files directly
    for file in jupiter_jrm09.dat jupiter_jrm33.dat mercury_thebault2018.dat mercury_wardinski2019.dat neptune_holme1996.dat saturn_cassini11.dat uranus_herbert2009.dat; do
    	curl -s "${base_url}/${file}" -o "data/${file}" 2>/dev/null && echo "Downloaded ${file}" || true
    done

libinternalfield:
    git clone https://github.com/mattkjames7/libinternalfield lib/libinternalfield
    cp -r lib/libinternalfield/data/coeffs data/coeffs