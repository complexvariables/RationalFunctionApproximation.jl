module RFAPythonCallExt

using RationalFunctionApproximation, PythonCall
const RFA = RationalFunctionApproximation

function RFA.approximate(f::Py, args...; kw...)
    try
        one = pyconvert(Vector, f)
        a  = [pyconvert(Vector, x) for x in args]
        RFA.approximate(one, a...; kw...)
    catch
        fun = x -> pyconvert(Number, pycall(f, x))
        RFA.approximate(fun, args...; kw...)
    end
end

end  # module
