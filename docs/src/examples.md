# Examples

## Computing

### Working with Units

By default, positions are in planetary radii

```@example example
using PlanetaryMagneticFields
using Unitful
model = JRM09()
r_km = 107238.0u"km"  # 1.5 RJ in km
θ, φ = π/4, 0.0
@assert model(r_km, θ, φ) == model(1.5, θ, φ)
model(r_km, θ, φ)
```

### Computing on multiple positions

```@example example
using PlanetaryMagneticFields

model = JRM09()

# Evaluate at different positions
positions = [
    (1.5, 0.0, 0.0),      # North pole, 1.5 RJ
    (2.0, π/2, 0.0),      # Equator, 2 RJ
    (1.0, π/4, π/2),      # 45° colatitude, 90° longitude
]
B = model.(positions)
```

### Computing on a Grid

You can compute field values over a latitude-longitude grid:

```@example example
# Compute radial field at 1.5 planetary radii
r = 1.5
Br = fieldmap(model, r, nlat=90, nlon=180; idx = 1)
# Field statistics
println("Br range: $(minimum(Br)) to $(maximum(Br)) nT")
```

### Using Cartesian Coordinates

```@example example
# Cartesian input and output
model = JRM09(in=:cartesian, out=:cartesian)
x, y, z = 1.0, 0.5, 0.5  # In Jupiter radii
B = model(x, y, z)

# Cartesian input, spherical output
model(x, y, z; out=:spherical)
```

## Visualization

```@example plotting
using CairoMakie, GeoMakie
using PlanetaryMagneticFields
```

### Different Projections

Various map projections are supported via PROJ strings:

```@example plotting
model = load_model(:earth, "igrf2020")

let fig = Figure(size=(1200, 400))
    sf1 = plot_fieldmap(fig[1, 1], model; axis=(; title="Mollweide"), dest="+proj=moll")
    sf2 = plot_fieldmap(fig[1, 2], model; axis=(; title="Robinson"), dest="+proj=robin")
    sf3 = plot_fieldmap(fig[1, 3], model; axis=(; title="Equal Earth"), dest="+proj=eqearth")
fig
end
```

### Plotting All Planets

Create a comprehensive figure showing magnetic field maps for all supported planets:

```@example plotting
plot_models(; r=0.9, figure = (; size=(1200, 400)), position = Makie.Right())
```

This displays the radial magnetic field (Bᵣ) at the surface for:
- **Mercury** (Anderson 2012 model)
- **Earth** (IGRF 2020 model)
- **Mars** (Langlais 2019 crustal field model)
- **Jupiter** (JRM09 model)
- **Saturn** (Cassini 11 model)
- **Uranus** (AH5 model)
- **Neptune** (GSFC O8 model)
- **Ganymede** (Kivelson 2002 model)

### Custom Multi-Planet Figures

You can create custom layouts with specific planets and models:

```@example plotting
# Giant planets comparison
planets = [
    (:Jupiter, "JRM09"),
    (:Saturn, "Cassini11"),
    (:Uranus, "AH5"),
    (:Neptune, "GSFCO8"),
]
models = map(planets) do (planet, model_name)
    load_model(planet, model_name)
end
```

```@example plotting
plot_models(; r=0.9, models = models[1:2])
```

```@example plotting
let fig = Figure()
    fpos = [(1, 1), (1, 2), (2, 1), (2,2)]
    for (idx, (planet, model_name)) in enumerate(planets)
        ax = GeoAxis(fig[fpos[idx]...];  title="$(planet) => $(model_name)")
        sf = plot_fieldmap!(ax, models[idx]; r=1.0)
    end
    Label(fig[0, :], "Giant Planet Magnetic Fields")
    fig
end
```

### Different Field Components

Plot different components of the magnetic field:

```@example plotting
using LinearAlgebra: norm

model = load_model(:jupiter, "jrm09")
let fig = Figure(; size=(600, 300))
    # Radial component
    plot_fieldmap(fig[1, 1], model; idx=1)
    # Field magnitude
    plot_fieldmap(fig[1, 2], model; idx=norm, colormap=:plasma)
    fig
end
```