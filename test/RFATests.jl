module RFATests

using RationalFunctionApproximation, ReTest, ComplexRegions, DoubleFloats, IntervalSets, Logging
const RFA = RationalFunctionApproximation

# Equivalent to `isapprox(f.(z), r.(z), norm=u->maximum(abs, u); kw...)`, written out
# because Julia 1.13.0-rc3's `isapprox` for arrays ignores the `norm` keyword
# (JuliaLang/LinearAlgebra.jl#1675).
function pass(f, r, z; atol::Real=0, rtol::Real=-1)
    fz, rz = f.(z), r.(z)
    if rtol < 0    # same default as `isapprox`
        rtol = atol > 0 ? 0 : sqrt(eps(float(real(promote_type(eltype(fz), eltype(rz))))))
    end
    nrm(u) = maximum(abs, u)
    return nrm(fz - rz) <= max(atol, rtol * max(nrm(fz), nrm(rz)))
end
logger = Logging.SimpleLogger(stderr, Logging.Error)
global_logger(logger)

include("aaa.jl")
include("discrete.jl")
include("discretized-path.jl")
include("real_interval.jl")
include("imag_interval.jl")
include("circle.jl")
include("custom.jl")
include("operations.jl")
include("parfrac.jl")
include("deprecated.jl")

end
