```@meta
CurrentModule = PlanetaryMagneticFields
```

# PlanetaryMagneticFields.jl

A Julia package for planetary magnetic field modeling.

## Overview

PlanetaryMagneticFields.jl provides a unified framework for working with magnetic field models of planets in our solar system. It implements spherical harmonic expansions with Schmidt semi-normalization, following standards used in geomagnetism and planetary science.

### Features

- **Multi-planetary support**: Currently Jupiter, with Earth, Saturn, and others coming soon
- **Multiple models**: JRM09, JRM33 for Jupiter
- **Flexible evaluation**: Support for both spherical and Cartesian coordinates

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
model_cart = JRM33(max_degree=13, in=:cartesian, out=:cartesian)
B = model_cart(1.0, 0.0, 0.5)  # Returns [B_x, B_y, B_z] in nT

# Keyword arguments take precedence over constructor arguments
B_sph = model_cart(1.0, 0.0, 0.5; out=:spherical)  # Returns [B_r, B_θ, B_φ]
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

## Examples

### Basic Field Evaluation

```@example example
using PlanetaryMagneticFields

# Load model
model = JRM09()

# Evaluate at different positions
positions = [
    (1.5, 0.0, 0.0),      # North pole, 1.5 RJ
    (2.0, π/2, 0.0),      # Equator, 2 RJ
    (1.0, π/4, π/2),      # 45° colatitude, 90° longitude
]
B = model.(positions)
```

### Using Cartesian Coordinates

```julia
# Cartesian input and output
model = JRM09(in=:cartesian, out=:cartesian)
x, y, z = 1.0, 0.5, 0.5  # In Jupiter radii
B = model(x, y, z)

# Cartesian input, spherical output
B_sph = model(x, y, z; out=:spherical)
```

### Working with Units

```@example example
# By default, positions are in planetary radii
using Unitful
model = JRM09()
r_km = 107238.0u"km"  # 1.5 RJ in km
θ, φ = π/4, 0.0
@assert model(r_km, θ, φ) == model(1.5, θ, φ)
model(r_km, θ, φ)
```

### Using Different Models

```@example example
# Load JRM09 (degree 10)
jrm09 = JRM09()

# Load JRM33 truncated to degree 13 (recommended)
jrm33 = JRM33(max_degree=13)

# Load full JRM33 (degree 30)
jrm33_full = JRM33()

# Compare models
r, θ, φ = 1.5, π/4, 0.0
B09 = jrm09(r, θ, φ)
B33 = jrm33(r, θ, φ)
```

## Available Models

List all available models

```@repl example
available_models()
```

### Jupiter

```@repl example
available_models(:jupiter)
```

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

```@repl comparison
using PythonCall
using PlanetaryMagneticFields
using Chairmarks
@py import JupiterMag as jm

r, θ, φ = 1.5, π/4, 0.0
jm.Internal.Config(Model="jrm33", CartesianIn=false, CartesianOut=false)
Br, Bt, Bp = jm.Internal.Field(r, θ, φ);
B_py = pyconvert(Vector{Float64}, [Br[0], Bt[0], Bp[0]])
model = JRM33(max_degree=13)
B = model(r, θ, φ)
@assert B_py ≈ B
@b model(r, θ, φ), jm.Internal.Field(r, θ, φ)
```

## API Reference

```@index
```

```@autodocs
Modules = [PlanetaryMagneticFields]
```
