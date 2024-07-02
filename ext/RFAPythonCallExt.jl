module RFAPythonCallExt

using RationalFunctionApproximation, PythonCall
const RFA = RationalFunctionApproximation

function RFA.aaa(f::Py, args...; kw...)
    RFA.aaa(x -> pyconvert(Number, pycall(f,x)), args...; kw...)
end

function RFA.approximate(f::Py, args...; kw...)
    RFA.approximate(x -> pyconvert(Number, pycall(f,x)), args...; kw...)
end

end  # module
