using Test
using MagneticModels
using LinearAlgebra

@testset "Model loading" begin
    # Test JRM09 loading
    jrm09 = MagneticModels.load_model("jrm09")
    @test jrm09.obj.radius == 71492.0
    @test uppercase(jrm09.name) == "JRM09"
    @test degree(jrm09) == 20
    @test order(jrm09) == 20

    # Test JRM33 loading with truncation
    jrm33_13 = MagneticModels.load_model("jrm33", max_degree = 13)
    @test degree(jrm33_13) == 13
    @test order(jrm33_13) == 13

    # Test full JRM33
    jrm33_full = MagneticModels.load_model("jrm33")
    @test degree(jrm33_full) == 30
end

@testset "Coefficient access" begin
    jrm09 = MagneticModels.load_model("jrm09")

    # Access dipole coefficient
    g10 = jrm09.coeffs[1, 0].g
    @test g10 != 0.0

    # https://ankitbarik.github.io/planetMagFields/examples.html
    @test jrm09.coeffs[2, 0].g == 11670.4
    @test jrm09.coeffs[4, 2].h == 27811.2
end

@testset "Available models" begin
    models = MagneticModels.available_models(:jupiter)
    @test "jrm09" in models
    @test "jrm33" in models
    @test length(models) >= 2
    @test all(m -> isa(m, String), models)
end


@testset "Model info" begin
    info = MagneticModels.jupiter_model_info("jrm09")
    @test info["name"] == "JRM09"
    @test haskey(info, "description")
    @test haskey(info, "reference")
    @test haskey(info, "doi")
    @test info["degree"] == 10

    info33 = MagneticModels.jupiter_model_info("jrm33")
    @test info33["degree"] == 30
    @test info33["recommended_degree"] == 13
end

@testset "JRM33 vs JRM09 comparison" begin
    jrm09 = MagneticModels.load_model("jrm09")
    jrm33 = MagneticModels.load_model("jrm33", max_degree = 10)

    # At the same truncation (degree 10), fields should be similar but not identical
    r, θ, φ = 2.0, π / 3, π / 4
    B_jrm09 = jrm09(r, θ, φ)
    B_jrm33 = jrm33(r, θ, φ)

    # Magnitudes should be similar (within ~20%)
    @test norm(B_jrm09) ≈ norm(B_jrm33) rtol = 0.01
end
