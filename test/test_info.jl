@testset "model_info" begin
    # From planet and model name
    info = model_info(:jupiter, "jrm09")
    @test isa(info, Dict{String, Any})
    @test haskey(info, "name")
    @test haskey(info, "description")
    @test haskey(info, "reference")
    @test info["name"] == "JRM09"

    # From loaded model (access underlying SphericalHarmonicModel)
    model = load_model(:jupiter, "jrm09")
    info2 = model_info(model.model)
    @test info2["name"] == info["name"]
end
