# MagneticModels.jl

[![Build Status](https://github.com/JuliaSpacePhysics/MagneticModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/MagneticModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/MagneticModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/MagneticModels.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

**A unified framework for planetary magnetic field modeling in Julia.**

MagneticModels.jl provides easy access to spherical harmonic models of planetary magnetic fields, with a focus on clean API design, performance, and extensibility.

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("MagneticModels")`

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg?logo=julia)](https://JuliaSpacePhysics.github.io/MagneticModels.jl/dev/)

## Features & Roadmap

- 🪐 **Multi-planetary support**: General type system and framework for magnetic field models
  - [x] Available astronomical objects: Jupiter, Earth, Saturn, Mercury, Mars, Ganymede
  - [ ] Model discovery API (available_models())
  - [ ] Model metadata system
- Testing
  - [x] Jupiter (JRM09, JRM33)
  - [ ] Earth models (IGRF, WMM)
  - [ ] Saturn models (Cassini-derived)
  - [ ] Mercury, Mars models
- 📐 **Flexible coordinates**: Spherical and Cartesian coordinate systems
- [ ] Model composition
- [ ] Time-dependent coefficients (secular variation)
- [ ] Performance optimizations (caching, pre-allocation, [SHTns.jl](https://github.com/fgerick/SHTns.jl))
  - [ ] Batch/vectorized evaluation
- [ ] Visualization extensions

## Quick Start

```julia
using MagneticModels

# Load a Jupiter magnetic field model by unique name
model = load_model(:JRM33; max_degree=13)
# Or use the convenience accessor
model = JRM33(max_degree=13)

# Evaluate the field at a position (1.5 Jupiter radii, 45° colatitude, 0° longitude)
r, θ, φ = 1.5, π/4, 0.0
B = model(r, θ, φ)  # Returns [B_r, B_θ, B_φ] in nanoTesla
```

## Elsewhere

- [MOP community code](https://lasp.colorado.edu/mop/missions/juno/community-code/): community code that may be useful to the wider Jovian community
- [planetMagFields](https://github.com/AnkitBarik/planetMagFields): access and analyze information about magnetic fields of planets in our solar system and visualize them in both 2D and 3D
- [JupiterMag](https://github.com/mattkjames7/JupiterMag): Python wrapper for a collection of Jovian magnetic field models written in C++ (see libjupitermag).
- [PSH](https://github.com/rjwilson-LASP/PSH): Planetary Spherical Harmonics community code
- [TsyganenkoModels.jl](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl): Tsyganenko's models for Earth's magnetosphere
- [Saturn-Mag-Model](https://github.com/NASA-Planetary-Science/Saturn-Mag-Model): FORTRAN source code for a Saturnian magnetospheric empirical magnetic field model derived from Cassini magnetometer data

- [SHTns](https://bitbucket.org/nschaeff/shtns/src/master/) & [SHTns.jl](https://github.com/fgerick/SHTns.jl): A high performance library for Spherical Harmonic Transform written in C

### Data sources

- [libinternalfield](https://github.com/mattkjames7/libinternalfield): a C++ library for various internal magnetic field models