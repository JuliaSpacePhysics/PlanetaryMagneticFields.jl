module PlanetaryMagneticFields

using StaticArrays
using LinearAlgebra
using LazyArrays
using TOML
using Dates

# Include core functionality
include("types.jl")
include("data.jl")
include("coefficients.jl")
include("cotrans.jl")
include("spherical_harmonics.jl")

# Include public API
include("api.jl")

# Export API
# export MagneticFieldModel, InternalFieldModel, ExternalFieldModel, SphericalHarmonicModel

# Main user-facing functions
export load_model, available_models, model_info, MagneticModel, fieldmap

# Convenience model accessors
export JRM09,
    JRM33,
    IGRF

# Utility functions that users might want
export degree, order

# Plotting functions (implemented in MakieExt)
export plot_fieldmap, plot_fieldmap!, plot_models

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

end # module PlanetaryMagneticFields
