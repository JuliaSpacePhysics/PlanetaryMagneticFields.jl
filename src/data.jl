function _available_models(path)
    data_dir = pkgdir(@__MODULE__, "data/coeffs", string(path))
    files = readdir(data_dir)
    # Extract model names from .dat files
    models = [splitext(f)[1] for f in files if endswith(f, ".dat")]
    return sort(models)
end

function available_models(body)
    bsym = Symbol(lowercase(string(body)))
    bsym in Celestial_bodies || error("Unsupported planet: $body. Currently supported: $Celestial_bodies")

    models = _available_models(bsym)
    if bsym == :earth
        append!(models, string.(TSYGANENKO_MODELS))
    end
    return sort(unique(models))
end

available_models() = reduce(vcat, available_models(body) for body in Celestial_bodies)
