const ComplexPath = ComplexRegions.AbstractPath
const ComplexCurve = ComplexRegions.AbstractCurve
const ComplexClosedPath = ComplexRegions.AbstractClosedPath
const ComplexClosedCurve = ComplexRegions.AbstractClosedCurve
const ComplexSCRegion = ComplexRegions.AbstractSimplyConnectedRegion
const RealComplex{T} = Union{T, ComplexValues.AnyComplex{T}}
const VectorRealComplex{T} = Union{Vector{T}, Vector{Complex{T}}}
const MatrixRealComplex{T} = Union{Matrix{T}, Matrix{Complex{T}}}
const VectorVectorRealComplex{T} = Union{Vector{Vector{T}},Vector{Vector{Complex{T}}}}

const unit_interval = Segment(-1.,1.)
const unit_circle = Circle(0., 1.)
const unit_disk = disk(0., 1.)
const Domain = Union{ComplexCurve, ComplexPath, ComplexSCRegion, AbstractVector}

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
