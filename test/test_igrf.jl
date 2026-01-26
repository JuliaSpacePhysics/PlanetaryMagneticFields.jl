using GeoCotrans
using Test
using Dates
using PlanetaryMagneticFields
import PlanetaryMagneticFields as PMF
model = PMF.IGRF()

@testset "Basic functionality" begin
    @test degree(model) == 13
    @test order(model) == 13
    @test length(PMF.epochs(model)) == 26  # 1900 to 2025 in 5-year increments
    @test PMF.epoch_range(model) == (Date(1900), Date(2025))
    @test PMF.is_time_varying(model) == true
end

@testset "Comparison with GeoCotrans at epoch points" begin
    test_years = Date.(1970:5:2015)
    test_times = [
        Date(2005, 7, 2),
        Date(2010, 6, 14),
        Date(2008, 3, 14),
        Date(2007, 9, 21),
    ]
    r = 1.0  # Earth radii
    θ = deg2rad(45)
    φ = deg2rad(45)
    for t in vcat(test_years, test_times)
        B_pmf = model(r, θ, φ, t)
        B_geo = collect(GeoCotrans.igrf_B(r * model.obj.radius, θ, φ, t))
        @test B_pmf ≈ B_geo rtol = 1.0e-8
    end
end

@testset "Coordinate system options" begin
    r, θ, φ = 1.0, deg2rad(45), deg2rad(45)
    t = Date(2020)

    B_sph = model(r, θ, φ, t)
    B_cart = model(r, θ, φ, t; out = :cartesian)
    using LinearAlgebra
    @test norm(B_sph) == norm(B_cart)
end

@testset "IGRF Static Model (backward compatibility)" begin
    model = load_model(:earth, "igrf2015")
    r = 1.0
    θ = deg2rad(45)
    φ = deg2rad(45)
    B = model(r, θ, φ)
    @test B ≈ collect(GeoCotrans.igrf_B(r * model.obj.radius, θ, φ, Date(2015)))
end
