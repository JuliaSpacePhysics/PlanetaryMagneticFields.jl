```@meta
CurrentModule = PlanetaryMagneticFields
```

# PlanetaryMagneticFields.jl

[![DOI](https://zenodo.org/badge/1140943453.svg)](https://doi.org/10.5281/zenodo.18428922)
[![version](https://juliahub.com/docs/General/PlanetaryMagneticFields/stable/version.svg)](https://juliahub.com/ui/Packages/General/PlanetaryMagneticFields)

A Julia package for planetary magnetic field modeling.

## Overview

PlanetaryMagneticFields.jl provides a unified framework for working with magnetic field models of planets in our solar system. It implements spherical harmonic expansions with Schmidt semi-normalization, following standards used in geomagnetism and planetary science.

### Features

- **Multi-planetary support**: Mercury, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, and Ganymede
- **Multiple models**: 80+ models available across all planets
- **Flexible evaluation**: Support for both spherical and Cartesian coordinates
- **Visualization**: Plot magnetic field maps using GeoMakie (Mollweide, Hammer projections)

## Installation

```julia
using Pkg
Pkg.add("PlanetaryMagneticFields")
```

## Quick Start

```@example quick_start
using PlanetaryMagneticFields

# Load a Jupiter magnetic field model by unique name
model = load_model(:JRM33; max_degree=13)

# Or use the convenience accessor
model = JRM33(max_degree=13)

# Evaluate the field at a position (in planetary radii)
# Position: 1.5 RJ, 45° colatitude, 0° longitude
r, θ, φ = 1.5, π/4, 0.0
B = model(r, θ, φ)  # Returns [B_r, B_θ, B_φ] in nT

# Use Cartesian coordinates
B = model(1.0, 0.0, 0.5; in=:cartesian)  # Returns [B_x, B_y, B_z] in nT

# Keyword arguments take precedence over constructor arguments
B_sph = model(1.0, 0.0, 0.5; in=:cartesian, out=:spherical)  # Returns [B_r, B_θ, B_φ]
```

## Coordinate Systems

PlanetaryMagneticFields.jl currently supports two coordinate systems:

### Spherical Coordinates

- `r`: Radial distance from planet center (in planetary radii or km)
- `θ`: Colatitude in radians [0, π] (0 at north pole, π at south pole)
- `φ`: East longitude in radians [0, 2π]

### Cartesian Coordinates

- `x`, `y`, `z`: Right-handed system with z-axis toward north pole
- x-axis points to (θ=π/2, φ=0)
- y-axis points to (θ=π/2, φ=π/2)

## Visualization

PlanetaryMagneticFields.jl provides visualization capabilities through GeoMakie. To use plotting functions, you need to load CairoMakie (or GLMakie) and GeoMakie:

```@example plotting
using CairoMakie, GeoMakie
using PlanetaryMagneticFields
```

Here is an example of how to plot a magnetic field map for a single planet:

```@example plotting
model = JRM09() # Load Jupiter's JRM09 model

# Create a field map with Mollweide projection
plot_fieldmap(model;
    r=1.0,
    axis=(; title="Jupiter Radial Magnetic Field (JRM09)"),
    dest="+proj=moll"
)
```

See [visualization examples](./examples.md#visualization) for more details.

## Available Models

List all available models:

```@example quick_start
available_models()
```

### Loading Models for Any Planet

```@example quick_start
# Load a model by planet and model name
earth_model = load_model(:earth, "igrf2020")
saturn_model = load_model(:saturn, "cassini11")
neptune_model = load_model(:neptune, "gsfco8")

# Evaluate at the surface
earth_model(1.0, π/4, 0.0)
```

### Mercury

```@example quick_start
available_models(:Mercury)
```

### Earth

```@example quick_start
available_models(:Earth)
```

The IGRF (International Geomagnetic Reference Field) models span from 1900 to 2025.

### Mars

```@example quick_start
available_models(:Mars)
```

Mars has crustal magnetic fields (no global dynamo).

### Jupiter

```@example quick_start
available_models(:Jupiter)
```

```@repl quick_start
model_info(:jrm33)
model_info("jrm09")
```

### Saturn

```@example quick_start
available_models(:Saturn)
```

### Uranus

```@example quick_start
available_models(:uranus)
```

### Neptune

```@example quick_start
available_models(:neptune)
```

### Ganymede

```@example quick_start
available_models(:ganymede)
```

## Physical Parameters

The package uses the following mean radii for unit conversions:

| Body | Mean Radius (km) |
|------|------------------|
| Mercury | 2,439.7 |
| Earth | 6,371.2 |
| Mars | 3,389.5 |
| Jupiter | 71,492.0 |
| Saturn | 60,268.0 |
| Uranus | 25,559.0 |
| Neptune | 24,764.0 |
| Ganymede | 2,634.1 |

## Mathematical Background

The magnetic scalar potential is expanded in spherical harmonics:

```math
V(r,θ,φ) = a \sum_{n=1}^{N} \sum_{m=0}^{n} \left(\frac{a}{r}\right)^{n+1}
           [g_n^m \cos(mφ) + h_n^m \sin(mφ)] P_n^m(\cos θ)
```

where:
- `a` is the planetary radius
- `g_n^m, h_n^m` are the Gauss coefficients (Schmidt semi-normalized)
- `P_n^m` are the associated Legendre polynomials (Schmidt semi-normalized)
- `(r,θ,φ)` are spherical coordinates

The magnetic field **B** = -∇V is computed from derivatives of the potential.

## References

1. Connerney, J. E. P., et al. (2018). A new model of Jupiter's magnetic field from Juno's first nine orbits. *Geophysical Research Letters*, 45(6), 2590-2596. doi:10.1002/2018GL077312

2. Connerney, J. E. P., et al. (2022). A new model of Jupiter's magnetic field at the completion of Juno's Prime Mission. *Journal of Geophysical Research: Planets*, 127(2), e2021JE007055. doi:10.1029/2021JE007055

3. Winch, D. E., et al. (2005). Geomagnetism and Schmidt quasi-normalization. *Geophysical Journal International*, 160(2), 487-504.

## Comparison
 
```@example comparison
using PythonCall
using PlanetaryMagneticFields
using Chairmarks
@py import JupiterMag as jm

r, θ, φ = 1.5, π/4, 0.0
jm.Internal.Config(Model="jrm33", CartesianIn=false, CartesianOut=false)
```

```@repl comparison
Br, Bt, Bp = jm.Internal.Field(r, θ, φ);
B_py = pyconvert(Vector{Float64}, [Br[0], Bt[0], Bp[0]])
model = JRM33(max_degree=13)
B = model(r, θ, φ)
@assert B_py ≈ B
@b $model($r, $θ, $φ), jm.Internal.Field($r, $θ, $φ)
```

## API Reference

```@index
```

```@autodocs
Modules = [PlanetaryMagneticFields]
```
