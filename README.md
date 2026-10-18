# PlanetaryMagneticFields.jl

[![DOI](https://zenodo.org/badge/1140943453.svg)](https://doi.org/10.5281/zenodo.18428922)

**A unified framework for planetary magnetic field modeling in Julia.**

PlanetaryMagneticFields.jl provides easy access to spherical harmonic models of planetary magnetic fields.

## Quick Start

```julia
using Pkg; Pkg.add("PlanetaryMagneticFields")
using PlanetaryMagneticFields

models = available_models()
# Load a Jupiter magnetic field model by unique name
model = load_model(:JRM33; max_degree=13)
# Or use the convenience accessor
model = JRM33(max_degree=13)

# Evaluate the field at a position (1.5 Jupiter radii, 45° colatitude, 0° longitude)
r, θ, φ = 1.5, π/4, 0.0
B = model(r, θ, φ)  # Returns [B_r, B_θ, B_φ] in nanoTesla

# Time-Varying Models (i.e. IGRF for Earth)
using Dates

model = IGRF()
model(r, θ, φ, DateTime(2020))
```

## Features

- 🪐 **Multi-planetary support**: General type system and framework for magnetic field models
  - [x] Available astronomical objects: Jupiter, Earth, Saturn, Mercury, Mars, Ganymede
    - [ ] Jupiter current sheet model

## Elsewhere

- [MOP community code](https://lasp.colorado.edu/mop/missions/juno/community-code/): community code that may be useful to the wider Jovian community
- [planetMagFields](https://github.com/AnkitBarik/planetMagFields): access and analyze information about magnetic fields of planets in our solar system and visualize them in both 2D and 3D
- [JupiterMag](https://github.com/mattkjames7/JupiterMag): Python wrapper for a collection of Jovian magnetic field models written in C++ (see libjupitermag).
- [PSH](https://github.com/rjwilson-LASP/PSH): Planetary Spherical Harmonics community code
- [TsyganenkoModels.jl](https://github.com/JuliaSpacePhysics/TsyganenkoModels.jl): Tsyganenko's models for Earth's magnetosphere
- [Saturn-Mag-Model](https://github.com/NASA-Planetary-Science/Saturn-Mag-Model): FORTRAN source code for a Saturnian magnetospheric empirical magnetic field model derived from Cassini magnetometer data

- [POT3D](https://github.com/predsci/POT3D): High Performance Potential Field Solver to approximate the solar coronal magnetic field using observed photospheric magnetic fields as a boundary condition
- [SHTns](https://bitbucket.org/nschaeff/shtns/src/master/) & [SHTns.jl](https://github.com/fgerick/SHTns.jl): A high performance library for Spherical Harmonic Transform written in C

### Data sources

- [libinternalfield](https://github.com/mattkjames7/libinternalfield): a C++ library for various internal magnetic field models
- [planetMagFields](https://github.com/AnkitBarik/planetMagFields): Mercury `wardinski2019` ([Wardinski et al. 2019](https://doi.org/10.1029/2018JE005835)) and Saturn `cassini11+` ([Cao et al. 2020](https://doi.org/10.1016/j.icarus.2019.113541))