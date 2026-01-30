module PlanetaryMagneticFieldsMakieExt

using PlanetaryMagneticFields
using GeoMakie
using GeoMakie.Makie
using LinearAlgebra: norm
using AxisKeys: axiskeys

using PlanetaryMagneticFields: _field_func
import PlanetaryMagneticFields: plot_fieldmap, plot_models


default_isvertical(pos) = pos isa Makie.Right || pos isa Makie.Left

function PlanetaryMagneticFields.plot_fieldmap!(
        ax, model;
        r = 1.0, nlat = 180, nlon = 360,
        idx = 1,
        colormap = :RdBu,
        shading = Makie.NoShading,
        kwargs...
    )

    field_map = fieldmap(model, r, nlat, nlon; idx)
    lons = axiskeys(field_map, 1)
    lats = axiskeys(field_map, 2)
    field = parent(field_map)
    # For diverging colormaps with Br, center at zero
    colorrange = if idx == 1
        maxabs = maximum(abs, field)
        (-maxabs, maxabs)
    else
        extrema(field)
    end
    return surface!(ax, lons, lats, field; colormap, colorrange, shading, kwargs...)
end

function PlanetaryMagneticFields.plot_fieldmap(model; r = 1.0, figure = (;), kwargs...)
    fig = Figure(; figure...)
    (ax, plot), _ = plot_fieldmap(fig[1, 1], model; r, kwargs...)
    return Makie.FigureAxisPlot(fig, ax, plot)
end

_label(idx) = if idx == 1
    "Bᵣ [nT]"
elseif idx == 2
    "Bθ [nT]"
elseif idx == 3
    "Bφ [nT]"
elseif idx == norm
    "|B| [nT]"
end

function PlanetaryMagneticFields.plot_fieldmap(
        fpos, model; r = 1.0,
        axis = (;),
        dest = "+proj=moll",
        idx = 1,
        position = Makie.Bottom(),
        vertical = default_isvertical(position),
        flipaxis = position isa Makie.Right,
        alignmode = Outside(),
        add_colorbar = true,
        kwargs...
    )

    ax = GeoAxis(fpos[1, 1]; dest, axis...)
    sf = plot_fieldmap!(ax, model; r = r, idx = idx, kwargs...)
    label = _label(idx)
    cb = add_colorbar ? Colorbar(fpos[1, 1, position], sf; label, vertical, flipaxis, alignmode) : nothing
    return Makie.AxisPlot(ax, sf), cb
end

function default_models()
    default_mappings = [
        :mercury => "anderson2012",
        :earth => "igrf2020",
        # :mars => "langlais2019", # too computationally expensive
        :jupiter => "jrm09",
        :saturn => "cassini11",
        :uranus => "ah5",
        :neptune => "gsfco8",
        # :ganymede => "kivelson2002a",
    ]

    return map(default_mappings) do (planet, model_name)
        load_model(planet, model_name)
    end
end

function plot_models(; r = 1.0, models = nothing, figure = (;), kwargs...)
    models = @something models default_models()
    # Create figure with 2 rows × 4 columns layout
    fig = Figure(; figure...)
    for (idx, model) in enumerate(models)
        row = (idx - 1) ÷ 3 + 1
        col = (idx - 1) % 3 + 1
        obj = model.obj
        title = "$(titlecase(string(obj.name))) ($(model.name))"
        axis = (; title)
        plot_fieldmap(GridLayout(fig[row, col]), model; r, axis, kwargs...)
    end
    return fig
end

end # module
