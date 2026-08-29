# Global registry mapping model names to (planet, model_string)
const MODEL_REGISTRY = Dict{Symbol,Tuple{Planet,String}}(
    :JRM09 => (Jupiter, "jrm09"),
    :JRM33 => (Jupiter, "jrm33"),
)
const Celestial_bodies = (:jupiter, :mercury, :earth, :mars, :ganymede, :neptune, :uranus, :saturn)

"""
    load_model(model; kwargs...)
    load_model(planet, model; kwargs...)

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
A callable object that can be invoked as `model(r, θ, φ; kwargs...)`

# Examples
```julia
# Load by unique model name
model = load_model(:JRM33; max_degree=13)
B = model(1.5, π/4, 0.0)  # Returns [B_r, B_θ, B_φ] in nT

# Load by planet and model name
model = load_model(:jupiter, "jrm09")
B = model(1.5, π/4, 0.0)

# With custom coordinate systems
model = load_model(:JRM33; max_degree=13, from=:cartesian, to=:cartesian)
B = model(1.0, 0.0, 0.5)  # Returns [B_x, B_y, B_z] in nT

# Override output coordinates at call time
B_sph = model(1.0, 0.0, 0.5; to=:spherical)  # Returns [B_r, B_θ, B_φ]
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

const _pkgdir = pkgdir(@__MODULE__)

function _load_model(p::Planet, model; max_degree=nothing, use_cache=false, kw...)
    data_path = "$_pkgdir/data/coeffs/$(p.name)/$(lowercase(model)).dat"
    coeffs = load_coefficients(data_path; use_cache)
    return SphericalHarmonicModel(model, coeffs, p; degree=max_degree, kw...)
end

"""
    fieldmap(model, r, nlat, nlon)
    fieldmap(model, r = 1.0; nlat = 180, nlon = 360)

Lazy magnetic-field grid at radial distance `r` [planetary radii]. Returns a
[`LazyFieldMap`](@ref) of size `(nlon, nlat)` whose elements `[Bᵣ, Bθ, Bφ]` are
computed on access. Project components by broadcasting; recover axes via
`field.lons` and `field.lats`.

# Example
```julia
model = load_model(:jupiter, "jrm09")
field = fieldmap(model; r=1.0)
Br = getindex.(field, 1)        # radial component
B  = norm.(field)               # magnitude
field.lons, field.lats          # axes in degrees
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

# IGRF model for testing with `GeoCotrans.IGRF`
function _IGRF(; kwargs...)
    tvc = load_igrf_epochs()
    return SphericalHarmonicModel("IGRF", tvc, Earth; frame=GEO())
end

for f in [:JRM09, :JRM33]
    @eval $f(; kwargs...) = load_model($(QuoteNode(f)); use_cache=true, kwargs...)
    @eval export $f
end



"""
    plot_fieldmap!(ax, model; r=1.0, nlat=180, nlon=360, kwargs...)

Plot a magnetic field map on an existing GeoAxis `ax`, with `nlat` and `nlon` determining the number of latitude and longitude points, respectively.

# Example
```julia
using CairoMakie, GeoMakie
using PlanetaryMagneticFields

fig = Figure()
ax = GeoAxis(fig[1,1]; dest="+proj=moll")
model = load_model(:jupiter, "jrm09")
sf = plot_fieldmap!(ax, model; r=1.0)
Colorbar(fig[1,2], sf; label="Br [nT]")
fig
```
"""
function plot_fieldmap! end

"""
    plot_fieldmap(model; r=1.0, dest="+proj=moll", kwargs...)

Create a complete figure with a magnetic field map.

# Arguments
- `model`: The magnetic field model
- `r::Real=1.0`: Radial distance in planetary radii

# Example
```julia
using CairoMakie, GeoMakie
using PlanetaryMagneticFields

model = load_model(:jupiter, "jrm09")
plot_fieldmap(model; r=1.0, title="Jupiter JRM09")
```
"""
function plot_fieldmap end

"""
    plot_models(; r=1.0, models=nothing, projection="+proj=moll", kwargs...)

Create a figure showing magnetic field maps of all available planets.

# Default Models
- Mercury: "anderson2012"
- Earth: "igrf2020"
- Mars: "langlais2019"
- Jupiter: "jrm09"
- Saturn: "cassini11"
- Uranus: "ah5"
- Neptune: "gsfco8"
- Ganymede: "kivelson2002a"
"""
function plot_models end