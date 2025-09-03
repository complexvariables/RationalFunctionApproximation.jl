module RFAForwardDiffExt

using RationalFunctionApproximation, ForwardDiff
const RFA = RationalFunctionApproximation

RFA.derivative(f::RFA.AbstractRationalFunction) = z -> derivative(f, z)

function RFA.derivative(f::RFA.AbstractRationalFunction, z::Number)
    ζ = [reim(complex(z))...]
    g = ForwardDiff.gradient(x -> real(f(complex(x...))), ζ)
    return conj(complex(g...))
end

end # module
