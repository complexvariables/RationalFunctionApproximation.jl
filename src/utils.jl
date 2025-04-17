const ComplexPath = ComplexRegions.AbstractPath
const ComplexCurve = ComplexRegions.AbstractCurve
const ComplexCurveOrPath = Union{ComplexCurve, ComplexPath}
const ComplexClosedPath = ComplexRegions.AbstractClosedPath
const ComplexClosedCurve = ComplexRegions.AbstractClosedCurve
const ComplexClosedCurveOrPath = Union{ComplexClosedCurve, ComplexClosedPath}
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

struct DiscretizedCurve{T,F}
    curve::T
    points::Matrix
    params::Matrix{F}
    next::Vector{Int}
end

function DiscretizedCurve(curve::ComplexCurveOrPath, s::AbstractVector; refinement=0, maxpoints=length(s))
    F = real_type(curve)
    params = Matrix{F}(undef, maxpoints, refinement+1)
    n = length(s)
    for i in 1:n-1
        δ = (s[i+1] - s[i]) / (refinement + 1)
        for j in 0:refinement
            params[i, j+1] = s[i] + j * δ
        end
    end
    params[n, 1] = s[end]
    points = Matrix{complex(F)}(undef, maxpoints, refinement+1)
    for i in 1:n-1, j in 1:refinement+1
        points[i, j] = curve(params[i, j])
    end
    points[n, 1] = curve(s[end])
    next = collect(2:n)
    sizehint!(next, maxpoints)
    return DiscretizedCurve{typeof(curve),F}(curve, points, params, next)
end

function DiscretizedCurve(curve::ComplexCurveOrPath, n::Integer; kwargs...)
    DiscretizedCurve(curve, range(0, length(curve), n); kwargs...)
end

function add_point!(d::DiscretizedCurve, idx)
    @assert length(idx) == 2
    (idx[2] == 1) && return
    n = length(d.next) + 1
    if n == size(d.points, 1)
        error("Cannot add more points to the discretization")
    end
    if idx[1] >= n
        error("Cannot add points beyond the last")
    end
    s_new = d.params[idx[1], idx[2]]
    ref = size(d.params, 2)
    # Replace row idx[1] with new values
    δ = (s_new - d.params[idx[1], 1]) / ref
    for j in 2:ref
        d.params[idx[1], j] = d.params[idx[1], 1] + (j - 1) * δ
        d.points[idx[1], j] = d.curve(d.params[idx[1], j])
    end
    # Add new row at the end
    succ = d.next[idx[1]]    # inherits the old successor
    δ = (d.params[succ, 1] - s_new) / ref
    for j in 1:ref
        d.params[n+1, j] = s_new + (j - 1) * δ
        d.points[n+1, j] = d.curve(d.params[n+1, j])
    end
    d.next[idx[1]] = n + 1    # new successor for the old point
    push!(d.next, succ)       # old succcessor for the new point
    return d.points[n+1, 1]
end

function gather(d::DiscretizedCurve, depth=0)
    n = length(d.next)
    points = Matrix{complex(real_type(d))}(undef, n, size(d.points, 2))
    params = Matrix{real_type(d)}(undef, n, size(d.params, 2))
    for i in 1:n
        points[i, :] = d.points[i, :]
        params[i, :] = d.params[i, :]
    end
    return points, params
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
