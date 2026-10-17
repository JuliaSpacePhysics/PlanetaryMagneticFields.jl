using GeoCotrans
using Test
using Dates
using PlanetaryMagneticFields

@testset "Earth IGRF files match GeoCotrans IGRF-14" begin
    r, θ, φ = 1.0, deg2rad(45), deg2rad(45)
    for y in 1965:5:2025
        @test load_model(:earth, "igrf$y")(r, θ, φ) ≈ GeoCotrans.igrf_B(r, θ, φ, Date(y)) rtol = 1.0e-12
    end
end
