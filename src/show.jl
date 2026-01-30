# Pretty printing
function Base.show(io::IO, m::SphericalHarmonicModel)
    print(io, "SphericalHarmonicModel(")
    print(io, "name=$(m.name), ")
    print(io, "obj=$(m.obj), ")
    print(io, "degree=$(degree(m)), ")
    print(io, "order=$(order(m))")
    if is_time_varying(m)
        t0, t1 = epoch_range(m)
        print(io, ", epochs=$(t0)-$(t1)")
    end
    print(io, ")")
    return
end

function Base.show(io::IO, coeffs::GaussCoefficients)
    return print(io, "GaussCoefficients(degree=$(degree(coeffs)), order=$(order(coeffs)))")
end
