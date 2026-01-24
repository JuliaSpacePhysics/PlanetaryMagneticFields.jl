```@meta
CurrentModule = MagneticModels
```

# MagneticModels.jl

A Julia package for planetary magnetic field modeling.

## Overview

MagneticModels.jl provides a unified framework for working with magnetic field models of planets in our solar system. It implements spherical harmonic expansions with Schmidt semi-normalization, following standards used in geomagnetism and planetary science.

### Features

- **Multi-planetary support**: Currently Jupiter, with Earth, Saturn, and others coming soon
- **Multiple models**: JRM09, JRM33 for Jupiter
- **Flexible evaluation**: Support for both spherical and Cartesian coordinates
- **Type-stable implementation**: Optimized for performance with StaticArrays
- **Simple API**: Easy to use for beginners, powerful for advanced users

## Installation

```julia
using Pkg
Pkg.add("MagneticModels")
```

## Quick Start

```julia
using MagneticModels

# Load a Jupiter magnetic field model
model = MagneticModel(:jupiter, "jrm33"; max_degree=13)

# Evaluate the field at a position (in planetary radii)
# Position: 1.5 RJ, 45° colatitude, 0° longitude

r, θ, φ = 8, π/2, 0.0
B = evaluate(model, [r, θ, φ])  # Returns [B_r, B_θ, B_φ] in nT

# Spherical coordinate example for scalar at 10 Rj, Colatitude on equator
# at East longitude of 38 degrees (converted to radians)
r, θ, φ = 10, π/2, deg2rad(38)
B = evaluate(model, [r, θ, φ])

# Get field magnitude
using LinearAlgebra
B_magnitude = norm(B)
println("Field magnitude: $B_magnitude nT")

# List available models
models = available_models(:jupiter)
println("Available models: ", join(models, ", "))
```

## Coordinate Systems

MagneticModels.jl supports two coordinate systems:

### Spherical Coordinates

- `r`: Radial distance from planet center (in planetary radii or km)
- `θ`: Colatitude in radians [0, π] (0 at north pole, π at south pole)
- `φ`: East longitude in radians [0, 2π]

### Cartesian Coordinates

- `x`, `y`, `z`: Right-handed system with z-axis toward north pole
- x-axis points to (θ=π/2, φ=0)
- y-axis points to (θ=π/2, φ=π/2)

## Examples

### Basic Field Evaluation

```julia
using MagneticModels

# Load model
model = MagneticModel(:jupiter, "jrm09")

# Evaluate at different positions
positions = [
    [1.5, 0.0, 0.0],      # North pole, 1.5 RJ
    [2.0, π/2, 0.0],      # Equator, 2 RJ
    [1.0, π/4, π/2],      # 45° colatitude, 90° longitude
]

for pos in positions
    B = evaluate(model, pos)
    println("Position: $pos → B = $B nT")
end
```

### Using Cartesian Coordinates

```julia
# Cartesian input and output
x, y, z = 1.0, 0.5, 0.5  # In Jupiter radii
B = evaluate(model, [x, y, z], coords=:cartesian, output_coords=:cartesian)
println("B = $B nT")

# Cartesian input, spherical output
B_sph = evaluate(model, [x, y, z], coords=:cartesian, output_coords=:spherical)
println("B (spherical) = $B_sph nT")
```

### Using Different Models

```julia
# Load JRM09 (degree 10)
jrm09 = MagneticModel(:jupiter, "jrm09")

# Load JRM33 truncated to degree 13 (recommended)
jrm33 = MagneticModel(:jupiter, "jrm33", max_degree=13)

# Load full JRM33 (degree 30)
jrm33_full = MagneticModel(:jupiter, "jrm33")

# Compare models
position = [1.5, π/4, 0.0]
B09 = evaluate(jrm09, position)
B33 = evaluate(jrm33, position)

println("JRM09 field magnitude: $(norm(B09)) nT")
println("JRM33 field magnitude: $(norm(B33)) nT")
```

### Model Information

```julia
# Get information about a model
info = model_info(:jupiter, "jrm09")
println("Model: ", info["name"])
println("Description: ", info["description"])
println("Reference: ", info["reference"])
println("DOI: ", info["doi"])
println("Degree: ", info["degree"])
```

### Working with Positions in Kilometers

```julia
# By default, positions are in planetary radii
# Use in_planetary_radii=false for positions in kilometers

r_km = 107238.0  # 1.5 RJ in km
θ, φ = π/4, 0.0

B = evaluate(model, [r_km, θ, φ], in_planetary_radii=false)
println("Field at $(r_km) km: $B nT")
```

## Available Models

### Jupiter

- **JRM09**: Juno Reference Model through Perijove 9 (degree 10)
  - Reference: Connerney et al. (2018), Geophysical Research Letters
  - DOI: 10.1002/2018GL077312

- **JRM33**: Juno Reference Model through Perijove 33 (degree 30, recommended truncation: 13)
  - Reference: Connerney et al. (2022), Journal of Geophysical Research: Planets
  - DOI: 10.1029/2021JE007055

## Physical Parameters

### Jupiter
- Mean radius: 71,492 km (1 RJ)
- Reference frame: System III (1965) coordinates
- Magnetic dipole moment: ~4.17 Gauss RJ³

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
@py import JupiterMag as jm

jm.Internal.Config(Model="jrm33",CartesianIn=false,CartesianOut=false)
Br,Bt,Bp = jm.Internal.Field(r, θ, φ)

# @py import planetmagfields: Planet
# p = Planet("jupiter")
# p.extrapolate([r])
# p.plot(r=r)
```

## API Reference

```@index
```

```@autodocs
Modules = [MagneticModels]
```
