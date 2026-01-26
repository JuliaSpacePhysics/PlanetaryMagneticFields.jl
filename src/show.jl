# Pretty printing
function Base.show(io::IO, model::SphericalHarmonicModel)
    print(io, "SphericalHarmonicModel(")
    print(io, "name=$(model.name), ")
    print(io, "degree=$(degree(model)), ")
    print(io, "order=$(order(model))")
    if is_time_varying(model)
        t0, t1 = epoch_range(model)
        print(io, ", epochs=$(t0)-$(t1)")
    end
    print(io, ")")
    return
end

function Base.show(io::IO, coeffs::GaussCoefficients)
    return print(io, "GaussCoefficients(degree=$(degree(coeffs)), order=$(order(coeffs)))")
end

# MagneticModel - check if underlying model is time-varying
function Base.show(io::IO, m::MagneticModel)
    print(io, "MagneticModel(")
    print(io, "name=$(m.name), ")
    print(io, "obj=$(m.obj), ")
    print(io, "input_coords=:$(m.in), ")
    print(io, "output_coords=:$(m.out), ")
    print(io, "model=$(m.model)")
    return print(io, ")")
end
