using CairoMakie
using GeoMakie
using PlanetaryMagneticFields
using Test

@testset "plot_fieldmap! on existing GeoAxis" begin
    model = JRM09()
    fig = Figure()
    ax = GeoAxis(fig[1, 1]; dest = "+proj=moll")
    sf = plot_fieldmap!(ax, model; r = 1.0, nlat = 18, nlon = 36)
    @test sf isa Makie.AbstractPlot
end

@testset "plot_fieldmap returns FigureAxisPlot" begin
    model = JRM09()
    fap = plot_fieldmap(model; r = 1.0, nlat = 18, nlon = 36)
    @test fap isa Makie.FigureAxisPlot
    @test fap.figure isa Figure
end

@testset "plot_fieldmap idx variants" begin
    using LinearAlgebra: norm
    model = JRM09()
    for idx in (1, 2, 3, norm)
        fap = plot_fieldmap(model; r = 1.0, nlat = 18, nlon = 36, idx)
        @test fap isa Makie.FigureAxisPlot
    end
end

@testset "plot_models with custom subset" begin
    models = [JRM09()]
    fig = plot_models(; r = 1.0, models, nlat = 18, nlon = 36)
    @test fig isa Figure
end
