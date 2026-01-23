using MagneticModels
using Test
using Aqua

@testset "MagneticModels.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(MagneticModels)
    end
    # Write your tests here.
end
