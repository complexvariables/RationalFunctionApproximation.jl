module RFAZygoteExt

using RationalFunctionApproximation, Zygote
const RFA = RationalFunctionApproximation

export derivative

function RFA.derivative(f::RFA.AbstractRationalFunction)
    return function(z::Number)
        ζ = complex(z)
        return something(conj(gradient(real ∘ f, ζ)[1]), zero(ζ))
    end
end

function RFA.derivative(f::RFA.AbstractRationalFunction, z::Number)
    ζ = complex(z)
    return something(conj(gradient(real ∘ f, ζ)[1]), zero(ζ))
end

end # module
