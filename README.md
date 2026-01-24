# MagneticModels.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSpacePhysics.github.io/MagneticModels.jl/dev/)
[![Build Status](https://github.com/JuliaSpacePhysics/MagneticModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/MagneticModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/MagneticModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/MagneticModels.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

**A unified framework for planetary magnetic field modeling in Julia.**

MagneticModels.jl provides easy access to spherical harmonic models of planetary magnetic fields, with a focus on clean API design, performance, and extensibility.

## Features

- 🪐 **Multi-planetary support**: Currently Jupiter (JRM09, JRM33), with more planets coming soon
- 🧮 **Spherical harmonic expansion**: Schmidt semi-normalized Legendre polynomials
- 📐 **Flexible coordinates**: Spherical and Cartesian coordinate systems
- ⚡ **High performance**: Type-stable implementation with StaticArrays
- 📚 **Well-documented**: Comprehensive documentation and examples
- 🧪 **Well-tested**: Extensive test suite with validation against reference implementations

## Installation

```julia
using Pkg
Pkg.add("MagneticModels")
```

## Quick Start

```julia
using MagneticModels

# Load a Jupiter magnetic field model
model = MagneticModel(:jupiter, "jrm33")

# Evaluate the field at a position (1.5 Jupiter radii, 45° colatitude, 0° longitude)
r, θ, φ = 1.5, π/4, 0.0
B = evaluate(model, [r, θ, φ])  # Returns [B_r, B_θ, B_φ] in nanoTesla

# Get field magnitude
using LinearAlgebra
B_magnitude = norm(B)
println("Magnetic field magnitude: $B_magnitude nT")
```

## Supported Models

### Jupiter
- **JRM09** (degree 10): Juno Reference Model through Perijove 9
  - Connerney et al. (2018), doi:10.1002/2018GL077312
- **JRM33** (degree 13-30): Juno Reference Model through Perijove 33
  - Connerney et al. (2022), doi:10.1029/2021JE007055

### Coming Soon
- Earth (IGRF-14, WMM2025)
- Saturn (Cassini models)
- Mercury, Mars, and other planets

## Examples

### Evaluate field in Cartesian coordinates

```julia
x, y, z = 1.0, 0.5, 0.5  # Position in Jupiter radii
B = evaluate(model, [x, y, z], coords=:cartesian, output_coords=:cartesian)
```

### Compare different models

```julia
jrm09 = MagneticModel(:jupiter, "jrm09")
jrm33 = MagneticModel(:jupiter, "jrm33", max_degree=13)

position = [1.5, π/4, 0.0]
B09 = evaluate(jrm09, position)
B33 = evaluate(jrm33, position)

println("JRM09: $(norm(B09)) nT")
println("JRM33: $(norm(B33)) nT")
```

### Get model information

```julia
info = model_info(:jupiter, "jrm09")
println(info["description"])  # "Juno Reference Model through Perijove 9"
println(info["reference"])    # Full citation
```

## Documentation

For detailed documentation, examples, and API reference, see the [documentation](https://JuliaSpacePhysics.github.io/MagneticModels.jl/dev/).

## Related Projects

## Elsewhere

- [MOP community code](https://lasp.colorado.edu/mop/missions/juno/community-code/): community code that may be useful to the wider Jovian community
- [planetMagFields](https://github.com/AnkitBarik/planetMagFields): access and analyze information about magnetic fields of planets in our solar system and visualize them in both 2D and 3D
- [PSH](https://github.com/rjwilson-LASP/PSH): Planetary Spherical Harmonics community code
- [TsyganenkoModels.jl](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl): Tsyganenko's models for Earth's magnetosphere
- [Saturn-Mag-Model](https://github.com/NASA-Planetary-Science/Saturn-Mag-Model): FORTRAN source code for a Saturnian magnetospheric empirical magnetic field model derived from Cassini magnetometer data

- [SHTns](https://bitbucket.org/nschaeff/shtns/src/master/) & [SHTns.jl](https://github.com/fgerick/SHTns.jl): A high performance library for Spherical Harmonic Transform written in C

### Data sources

- [libinternalfield](https://github.com/mattkjames7/libinternalfield): a C++ library for various internal magnetic field models