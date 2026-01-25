function _available_models(path)
    data_dir = pkgdir(@__MODULE__, "data/coeffs", string(path))
    files = readdir(data_dir)
    # Extract model names from .dat files
    models = [splitext(f)[1] for f in files if endswith(f, ".dat")]
    return sort(models)
end

function available_models(body)
    bsym = Symbol(lowercase(string(body)))
    return if bsym in Celestial_bodies
        _available_models(bsym)
    else
        error("Unsupported planet: $body. Currently supported: $Celestial_bodies")
    end
end

available_models() = reduce(vcat, available_models(body) for body in Celestial_bodies)
