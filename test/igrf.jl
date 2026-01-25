using GeoCotrans
using Test
using Dates

@testset "IGRF" begin
    model = load_model(:earth, "igrf2015")
    r = 1.0
    θ = deg2rad(45)
    φ = deg2rad(45)
    B = model(r, θ, φ)
    @test B ≈ collect(GeoCotrans.igrf_B(r * model.obj.radius, θ, φ, Date(2015)))
end
