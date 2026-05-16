using BenchmarkTools
using PlanetaryMagneticFields
using Dates

const SUITE = BenchmarkGroup()

# Model loading
SUITE["load"] = BenchmarkGroup()
SUITE["load"]["JRM09"]  = @benchmarkable JRM09()
SUITE["load"]["JRM33"]  = @benchmarkable JRM33()
SUITE["load"]["JRM33_d13"] = @benchmarkable JRM33(max_degree = 13)
SUITE["load"]["IGRF"]   = @benchmarkable IGRF()

# Single-point evaluation
SUITE["eval"] = BenchmarkGroup()
let jrm33 = JRM33(), jrm33_13 = JRM33(max_degree = 13), igrf = IGRF()
    r, θ, φ = 1.5, π / 4, 0.0
    SUITE["eval"]["JRM33"]     = @benchmarkable $jrm33($r, $θ, $φ)
    SUITE["eval"]["JRM33_d13"] = @benchmarkable $jrm33_13($r, $θ, $φ)
    SUITE["eval"]["IGRF"]      = @benchmarkable $igrf($r, $θ, $φ, $(Date(2015)))
end

# Field map over lat-lon grid
SUITE["fieldmap"] = BenchmarkGroup()
let jrm33_13 = JRM33(max_degree = 13)
    SUITE["fieldmap"]["JRM33_d13_45x90"]   = @benchmarkable fieldmap($jrm33_13; nlat = 45,  nlon = 90)
    SUITE["fieldmap"]["JRM33_d13_180x360"] = @benchmarkable fieldmap($jrm33_13; nlat = 180, nlon = 360)
end
