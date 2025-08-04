module RFATests

using RationalFunctionApproximation, ReTest, ComplexRegions, DoubleFloats, Logging
const RFA = RationalFunctionApproximation

pass(f, r, z; kw...) = isapprox(f.(z), r.(z), norm=u->maximum(abs, u); kw...)
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

end
