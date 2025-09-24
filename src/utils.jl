const ComplexPath{T} = ComplexRegions.AbstractPath{T}
const ComplexCurve{T} = ComplexRegions.AbstractCurve{T}
const ComplexCurveOrPath{T} = Union{ComplexCurve{T}, ComplexPath{T}}
const ComplexClosedPath{T} = ComplexRegions.AbstractClosedPath{T}
const ComplexClosedCurve{T} = ComplexRegions.AbstractClosedCurve{T}
const ComplexClosedCurveOrPath{T} = Union{ComplexClosedCurve{T}, ComplexClosedPath{T}}
const ComplexSCRegion = ComplexRegions.AbstractSimplyConnectedRegion
const RealComplex{T} = Union{T, ComplexValues.AnyComplex{T}}
const VectorRealComplex{T} = Union{Vector{T}, Vector{Complex{T}}}
const MatrixRealComplex{T} = Union{Matrix{T}, Matrix{Complex{T}}}
const VectorVectorRealComplex{T} = Union{Vector{Vector{T}},Vector{Vector{Complex{T}}}}

const unit_interval = Segment(-1.,1.)
const unit_circle = Circle(0., 1.)
const unit_disk = disk(0., 1.)
const Domain{T} = Union{ComplexCurve{T}, ComplexPath{T}, ComplexSCRegion{T}, AbstractVector{T}}

isclosed(p::ComplexCurve) = isa(p, ComplexClosedCurve)
isclosed(p::ComplexPath) = isa(p, ComplexClosedPath)

# Provide a fallback `dist` method for unknown (i.e., user-defined) curve types
function ComplexRegions.dist(z::Number, p::ComplexRegions.AbstractCurve)
    return minimum(abs(z - point(p, t)) for t in range(0, 1, length=300))
end

# Refinement in parameter space
function refine(t, N, include_ends=false)
    x = sort(t)
    Δx = diff(x)
    if include_ends
        d = eltype(x).((0:N) / (N+1))
    else
        d = eltype(x).((1:N) / (N+1))
    end
    return vec( x[1:end-1] .+ (d' .* Δx) )
end

# Refinement on a curve
refine(p::ComplexCurve, args...) = refine(Path(p), args...)
function refine(p::ComplexPath, t::AbstractVector, N::Integer=3, args...)
    tt = refine(t, N, args...)
    ττ = point(p, tt)
    if isreal(p)
        ττ = real(ττ)
    end
    return tt, ττ
end

function golden_max(f::Function, a::Number, b::Number; tol=sqrt(eps(float(b-a))))
    a, b, five = promote(float(a), float(b), 5)
    φ = (sqrt(five) - 1) / 2
    φʹ = 1 - φ
    x1 = a + φʹ * (b - a)
    x2 = a + φ * (b - a)
    f1 = f(x1)
    f2 = f(x2)
    while b - a > 2tol
        if f1 > f2
            b = x2
            x2, f2 = x1, f1
            x1 = a + φʹ * (b - a)
            f1 = f(x1)
        else
            a = x1
            x1, f1 = x2, f2
            x2 = a + φ * (b - a)
            f2 = f(x2)
        end
    end
    return (a + b) / 2
end
