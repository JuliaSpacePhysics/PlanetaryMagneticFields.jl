"""
Public API for PlanetaryMagneticFields.jl

This module provides the main user-facing functions for loading and evaluating
magnetic field models.
"""

# Global registry mapping model names to (planet, model_string)
const MODEL_REGISTRY = Dict{Symbol, Tuple{Symbol, String}}(
    :JRM09 => (:jupiter, "jrm09"),
    :JRM33 => (:jupiter, "jrm33"),
    :T89 => (:earth, "t89"),
    :T96 => (:earth, "t96"),
    :T01 => (:earth, "t01"),
    :TS04 => (:earth, "ts04"),
)
const Celestial_bodies = (:jupiter, :mercury, :earth, :mars, :ganymede, :neptune, :uranus, :saturn)
const DEFAULT_IGRF = "igrf2025"

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

load_model(p, model; kw...) = _load_model(planet(p), model; kw...)

function _latest_igrf_model()
    igrf_models = filter(name -> startswith(name, "igrf"), _available_models(:earth))
    isempty(igrf_models) && return DEFAULT_IGRF
    years = map(name -> parse(Int, replace(name, "igrf" => "")), igrf_models)
    latest_idx = findmax(years)[2]
    return igrf_models[latest_idx]
end

function _load_model(
        p::Planet,
        model;
        max_degree = nothing,
        drivers = nothing,
        igrf_model = nothing,
        in = :spherical,
        out = :spherical,
    )
    model_str = lowercase(String(model))
    model_sym = Symbol(model_str)

    if p.name == :earth && Base.in(model_sym, TSYGANENKO_MODELS)
        drv = isnothing(drivers) ? TsyganenkoDrivers() : drivers
        drv isa TsyganenkoDrivers || error("drivers must be a TsyganenkoDrivers struct")
        igrf_choice = isnothing(igrf_model) ? _latest_igrf_model() : igrf_model
        internal = _load_model(p, igrf_choice; max_degree = max_degree, in = :spherical, out = :spherical)
        external = tsyg_model(p, model_sym, drv)
        return compose_models((internal, external); in = in, out = out)
    end

    data_file = "$(model_str).dat"
    data_path = pkgdir(@__MODULE__, "data/coeffs", lowercase(string(p.name)), data_file)
    coeffs = load_coefficients(data_path, max_degree = max_degree)
    sh_model = SphericalHarmonicModel(uppercase(model_str), coeffs)
    return MagneticModel(sh_model, p; in = in, out = out)
end

"""
    fieldmap(model, r, nlat, nlon; idx = identity)
    fieldmap(model, r = 1.0; nlat = 180, nlon = 360, kw...)

Compute magnetic fields over a latitude-longitude grid at a given radial distance `r` [planetary radii] for model `model`.

# Arguments
- `nlat=180`: Number of latitude points
- `nlon=360`: Number of longitude points`
- `idx`: Field component selector / function that works on the output of `model(r, θ, φ)` (default: identity returns full vector)
  - `1`: Radial component (Br)
  - `2`: Colatitude component (Bθ)
  - `3`: Azimuthal component (Bφ)
  - `norm`: Field magnitude

# Returns
- `KeyedArray`: Magnetic field values at each grid point with axes:
  - `lon`: Longitude values [-180, 180] degrees
  - `lat`: Latitude values [-90, 90] degrees

# Example
```julia
model = load_model(:jupiter, "jrm09")
field_map = fieldmap(model; r=1.0)
# Access data: field_map[lon=0.0, lat=45.0]
# Get axes: axiskeys(field_map, 1) for longitudes, axiskeys(field_map, 2) for latitudes
```
"""
function fieldmap end

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
    model_info(model)

Get detailed information about a magnetic field model.

# Example
```julia
info = model_info("jrm09")
println(info["description"])
```
"""
function model_info end

function model_info(model)
    name = lowercase(String(model))
    file = pkgdir(@__MODULE__, "data/metadata/$(name).toml")
    @assert isfile(file)
    return TOML.parsefile(file)
end

for f in [:JRM09, :JRM33]
    @eval $f(; kwargs...) = load_model($(QuoteNode(f)); kwargs...)
    @eval export $f
end
