using MagneticModels
using Test
using Aqua
using MagneticModels: SphericalHarmonicModel, GaussCoefficients
using Unitful

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(MagneticModels)
end

@testset "quick test" begin
    model = JRM33(max_degree = 13)
    r, θ, φ = 8, π / 2, 0.0
    B = model(r, θ, φ)
    @test B ≈ [-250.03964154, 779.36280353, -48.0067748]
    @test model(8 * 71492.0u"km", θ, φ) ≈ [-250.03964154, 779.36280353, -48.0067748]
end

@testset "MagneticModels.jl" begin
    # Test type system
    include("test_types.jl")
end

@testset "Spherical Harmonics" begin
    include("test_spherical_harmonics.jl")
end

@testset "Public API" begin
    include("test_api.jl")
end

@testset "Jupiter Models" begin
    include("test_jupiter.jl")
end
