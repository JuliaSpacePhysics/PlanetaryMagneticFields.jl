using Test
using PlanetaryMagneticFields

@testset "GaussCoefficients construction and indexing" begin
    # Create simple test coefficients (degree 2, order 2)
    g = zeros(3, 3)  # Rows: degree 0,1,2; Cols: order 0,1,2
    h = zeros(3, 3)

    # Set some values
    g[2, 1] = 410244.7  # g_1^0
    g[2, 2] = -71498.3  # g_1^1
    h[2, 2] = 21330.5   # h_1^1

    coeffs = GaussCoefficients(g, h)

    @test degree(coeffs) == 2
    @test order(coeffs) == 2

    # Test indexing
    c = coeffs[1, 0]
    @test c.g == 410244.7
    @test c.h == 0.0

    c = coeffs[1, 1]
    @test c.g == -71498.3
    @test c.h == 21330.5
    # Test bounds checking
    @test_throws BoundsError coeffs[3, 0]  # Degree too high
end

@testset "GaussCoefficients validation" begin
    g = zeros(3, 3)
    h = zeros(3, 3)

    check = true
    # Valid construction
    @test isa(GaussCoefficients(g, h; check = check), GaussCoefficients)
    # Invalid: mismatched sizes
    @test_throws ErrorException GaussCoefficients(g, zeros(2, 2); check = check)
    # Invalid: order > degree (more columns than rows)
    @test_throws ErrorException GaussCoefficients(zeros(2, 3), zeros(2, 3); check = check)
    # Invalid: degree < 1 (matrix too small)
    @test_throws ErrorException GaussCoefficients(zeros(1, 1), zeros(1, 1); check = check)
end

@testset "Pretty printing" begin
    g = zeros(3, 3)
    h = zeros(3, 3)
    coeffs = GaussCoefficients(g, h)

    # Test that show methods work without errors
    io = IOBuffer()
    show(io, coeffs)
    @test occursin("GaussCoefficients", String(take!(io)))

    model = SphericalHarmonicModel("test", coeffs)
    show(io, model)
    output = String(take!(io))
    @test occursin("SphericalHarmonicModel", output)
    @test occursin("test", output)
end
