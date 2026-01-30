using Test
using PlanetaryMagneticFields
using LinearAlgebra

@testset "CON2020 Model Loading" begin
    model = CON2020()
    @test model.name == "CON2020"
    @test model.obj.name == :jupiter
    @test model.obj.radius == 71492.0

    # Test model parameters
    @test model.params.r0 == 7.8
    @test model.params.r1 == 51.4
    @test model.params.d == 3.6
    @test model.params.μ_i_half == 139.6
    @test model.params.I_rho == 16.7
    @test model.params.use_radial_current == true
end

@testset "CAN1981 Model Loading" begin
    model = CAN1981()
    @test model.params.use_radial_current == false
    @test model.params.r0 == 5.0
    @test model.params.r1 == 50.0
end

@testset "CON2020 Field Evaluation" begin
    model = CON2020()

    r, θ, φ = 10.0, π / 2, 0.0
    B = model(r, θ, φ)

    @test length(B) == 3
    @test all(isfinite.(B))
    @test norm(B) > 0

    # Test at different radii
    for r_test in [8.0, 15.0, 30.0, 50.0]
        B_test = model(r_test, π / 2, 0.0)
        @test all(isfinite.(B_test))
    end
end

@testset "CON2020 Symmetry Properties" begin
    model = CON2020(θ_tilt=0.0, φ_tilt=0.0)

    r, θ = 15.0, π / 2
    B1 = model(r, θ, 0.0)
    B2 = model(r, θ, π / 2)
    B3 = model(r, θ, π)

    # Br and Bθ should be same for all φ (axisymmetric with no tilt)
    @test B1[1] ≈ B2[1] rtol = 0.01
    @test B1[1] ≈ B3[1] rtol = 0.01
    @test B1[2] ≈ B2[2] rtol = 0.01

    # Test north-south antisymmetry in Br
    θ_north = π / 2 - 0.1
    θ_south = π / 2 + 0.1

    B_north = model(r, θ_north, 0.0)
    B_south = model(r, θ_south, 0.0)
    @test B_north[1] ≈ -B_south[1] rtol = 0.1
end

@testset "CON2020 Equation Types" begin
    r, θ, φ = 15.0, π / 2, 0.0

    model_analytical = CON2020(equation_type=:analytical)
    model_integral = CON2020(equation_type=:integral)
    model_hybrid = CON2020(equation_type=:hybrid)

    B_analytical = model_analytical(r, θ, φ)
    B_integral = model_integral(r, θ, φ)
    B_hybrid = model_hybrid(r, θ, φ)

    @test all(isfinite.(B_analytical))
    @test all(isfinite.(B_integral))
    @test all(isfinite.(B_hybrid))

    # Results should be similar
    @test norm(B_analytical) ≈ norm(B_integral) rtol = 0.3
    @test norm(B_hybrid) ≈ norm(B_integral) rtol = 0.3
end

@testset "CON2020 Custom Parameters" begin
    model_custom = CON2020(
        μ_i_half=100.0,
        r0=6.0,
        r1=40.0,
        d=3.0,
        θ_tilt=10.0,
        φ_tilt=160.0,
        I_rho=15.0
    )

    @test model_custom.params.μ_i_half == 100.0
    @test model_custom.params.r0 == 6.0
    @test model_custom.params.r1 == 40.0

    B = model_custom(12.0, π / 2, 0.0)
    @test all(isfinite.(B))
end

@testset "CON2020 Radial Current" begin
    r, θ, φ = 15.0, π / 2 - 0.1, 0.0

    model_with_ir = CON2020(use_radial_current=true)
    model_without_ir = CON2020(use_radial_current=false)

    B_with = model_with_ir(r, θ, φ)
    B_without = model_without_ir(r, θ, φ)

    # Fields should be different due to radial current
    @test B_with != B_without
    @test B_with[1] ≈ B_without[1] rtol = 0.3
end

@testset "CON2020 + JRM33 Combined Field" begin
    jrm33 = JRM33(max_degree=13)
    con2020 = CON2020()

    r, θ, φ = 15.0, π / 2, 0.0

    B_internal = jrm33(r, θ, φ)
    B_external = con2020(r, θ, φ)
    B_total = B_internal .+ B_external

    @test all(isfinite.(B_total))
    @test norm(B_external) > 0.1 * norm(B_internal)
end

@testset "CON2020 Edge Cases" begin
    model = CON2020()

    # Near the axis
    B_axis = model(10.0, 0.01, 0.0)
    @test all(isfinite.(B_axis))

    # At inner edge
    B_inner = model(7.8, π / 2, 0.0)
    @test all(isfinite.(B_inner))

    # At outer edge
    B_outer = model(51.4, π / 2, 0.0)
    @test all(isfinite.(B_outer))

    # Outside the disc
    B_outside = model(60.0, π / 2, 0.0)
    @test all(isfinite.(B_outside))

    # Inside inner edge
    B_inside = model(5.0, π / 2, 0.0)
    @test all(isfinite.(B_inside))
end
