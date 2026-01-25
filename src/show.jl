# Pretty printing
function Base.show(io::IO, model::SphericalHarmonicModel)
    print(io, "SphericalHarmonicModel(")
    print(io, "name=$(model.name), ")
    print(io, "degree=$(degree(model)), ")
    print(io, "order=$(order(model))")
    print(io, ")")
    return
end

function Base.show(io::IO, coeffs::GaussCoefficients)
    return print(io, "GaussCoefficients(degree=$(coeffs.degree), order=$(coeffs.order))")
end


function Base.show(io::IO, m::MagneticModel)
    print(io, "MagneticModel(")
    print(io, "name=$(m.name), ")
    print(io, "obj=$(m.obj), ")
    print(io, "degree=$(degree(m.model)), ")
    print(io, "input_coords=:$(m.in), ")
    print(io, "output_coords=:$(m.out)")
    return print(io, ")")
end
