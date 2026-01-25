"""
Public API for MagneticModels.jl

This module provides the main user-facing functions for loading and evaluating
magnetic field models.
"""

# Global registry mapping model names to (planet, model_string)
const MODEL_REGISTRY = Dict{Symbol, Tuple{Symbol, String}}(
    :JRM09 => (:jupiter, "jrm09"),
    :JRM33 => (:jupiter, "jrm33"),
)
const Celestial_bodies = (:jupiter, :mercury, :earth, :mars, :ganymede, :neptune, :uranus, :saturn)

"""
    load_model(model; kwargs...)::MagneticModel
    load_model(planet, model; kwargs...)::MagneticModel

Load and create a magnetic field model.

# Single-argument form (by unique model name)
Load a model by its unique name (e.g., `:JRM33`, `:JRM09`).

# Two-argument form (by planet and model)
Load a model by specifying both planet and model name.

# Arguments
- `model`: Unique model identifier (e.g., `:JRM33`, `:JRM09`)
- `planet`: Planet identifier (e.g., :jupiter, :saturn, :earth)
- `max_degree`: Maximum degree to use (truncates model if specified)
- `in`: Default input coordinate system (`:spherical` or `:cartesian`)
- `out`: Default output coordinate system (`:spherical` or `:cartesian`)

# Notes
Input positions are assumed to be in planetary radii by default. For physical units (e.g., km),
use Unitful.jl which is supported via a package extension.

# Returns
A callable `MagneticModel` object that can be invoked as `model(r, θ, φ; kwargs...)`

# Examples
```julia
# Load by unique model name
model = load_model(:JRM33; max_degree=13)
B = model(1.5, π/4, 0.0)  # Returns [B_r, B_θ, B_φ] in nT

# Load by planet and model name
model = load_model(:jupiter, "jrm09")
B = model(1.5, π/4, 0.0)

# With custom coordinate systems
model = load_model(:JRM33; max_degree=13, in=:cartesian, out=:cartesian)
B = model(1.0, 0.0, 0.5)  # Returns [B_x, B_y, B_z] in nT

# Override output coordinates at call time
B_sph = model(1.0, 0.0, 0.5; out=:spherical)  # Returns [B_r, B_θ, B_φ]
```
"""
function load_model(name; kw...)
    # Look up planet and model string from registry
    key = name isa String ? Symbol(uppercase(name)) : name
    if !haskey(MODEL_REGISTRY, key)
        error("Unknown model: $name. Available models: $(keys(MODEL_REGISTRY))")
    end
    planet, model_str = MODEL_REGISTRY[key]
    return load_model(planet, model_str; kw...)
end

function _load_model(p::Planet, model; max_degree = nothing, kw...)
    data_file = "$(lowercase(model)).dat"
    data_path = pkgdir(@__MODULE__, "data/coeffs", lowercase(string(p.name)), data_file)
    coeffs = load_coefficients(data_path, max_degree = max_degree)
    sh_model = SphericalHarmonicModel(uppercase(model), coeffs)
    return MagneticModel(sh_model, p; kw...)
end

load_model(p, model; kw...) = _load_model(planet(p), model; kw...)

"""
    available_models([body])

List all available models, optionally specified by a given astronomical `body`.

# Example
```julia
available_models()
available_models(:jupiter)
```
"""
function available_models end

"""
    model_info(planet::Symbol, model::String)
    model_info(model::SphericalHarmonicModel)

Get detailed information about a magnetic field model.

# Arguments
- `planet::Symbol`: Planet identifier
- `model::String`: Model identifier
OR
- `model::SphericalHarmonicModel`: An already-loaded model

# Returns
- `Dict{String,Any}`: Model metadata

# Example
```julia
info = model_info(:jupiter, "jrm09")
println(info["description"])
println("DOI: ", info["doi"])

# Or from a loaded model
model = MagneticModel(:jupiter, "jrm09")
info = model_info(model)
```
"""
function model_info(planet::Symbol, model_name::String)
    if planet == :jupiter
        return jupiter_model_info(model_name)
    else
        error("Unsupported planet: $planet")
    end
end

function model_info(model::SphericalHarmonicModel)
    return model_info(model.planet, lowercase(model.name))
end

for f in [:JRM09, :JRM33]
    @eval $f(; kwargs...) = load_model($(QuoteNode(f)); kwargs...)
    @eval export $f
end
