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

struct DiscretizedPath{T,F}
    path::T
    points::Matrix
    params::Matrix{F}
    next::Vector{Int}
end

function DiscretizedPath(path::ComplexCurveOrPath, s::AbstractVector; refinement=0, maxpoints=length(s))
    F = real_type(path)
    params = Matrix{F}(undef, maxpoints, refinement+1)
    n = length(s)
    for i in 1:n-1
        δ = (s[i+1] - s[i]) / (refinement + 1)
        for j in 0:refinement
            params[i, j+1] = s[i] + j * δ
        end
    end
    points = Matrix{typeof(point(path, 0))}(undef, maxpoints, refinement+1)
    for i in 1:n-1, j in 1:refinement+1
        points[i, j] = path(params[i, j])
    end
    next = [collect(2:n-1); 0]
    sizehint!(next, maxpoints)
    return DiscretizedPath{typeof(path),F}(path, points, params, next)
end

function DiscretizedPath(path::ComplexCurveOrPath=ComplexRegions.Circle(0,1), n::Integer=0; kwargs...)
    DiscretizedPath(path, range(0, length(path), n); kwargs...)
end

function add_point!(d::DiscretizedPath, idx)
    @assert length(idx) == 2
    (idx[2] == 1) && return
    n = length(d.next) + 1
    if n == size(d.points, 1)
        error("Cannot add more points to the discretization")
    end
    # if d.next[idx[1]] == 0
    #     error("Cannot add points beyond the last")
    # end
    s_new = d.params[idx[1], idx[2]]
    ref = size(d.params, 2)
    # Replace row idx[1] with new values
    δ = (s_new - d.params[idx[1], 1]) / ref
    for j in 2:ref
        d.params[idx[1], j] = d.params[idx[1], 1] + (j - 1) * δ
        d.points[idx[1], j] = point(d.path, d.params[idx[1], j])
    end
    # Add new row at the end
    succ = d.next[idx[1]]    # inherits the old successor
    if succ == 0   # successor is the end of the path
        δ = (length(d.path) - s_new) / ref
    else
        δ = (d.params[succ, 1] - s_new) / ref
    end
    for j in 1:ref
        d.params[n, j] = s_new + (j - 1) * δ
        d.points[n, j] = point(d.path, d.params[n, j])
    end
    d.next[idx[1]] = n        # new successor for the old point
    push!(d.next, succ)       # old succcessor for the new point
    # return indexes into the new parts
    return [CartesianIndices((idx[1]:idx[1], 1:ref)); CartesianIndices((n:n, 1:ref))]
    # return view(d.params, newbies, :), view(d.points, newbies, :)
end

function Base.collect(d::DiscretizedPath, depth=0)
    n = length(d.next)
    if depth == :all
        depth = size(d.params, 2) - 1
    end
    params = similar(d.params, depth+1, n)
    points = similar(d.points, depth+1, n)
    k = idx = 1
    while idx != 0
        params[:, k] .= d.params[idx, 1:depth+1]
        points[:, k] .= d.points[idx, 1:depth+1]
        idx = d.next[idx]
        k += 1
    end
    s_fin = length(d.path)
    return [vec(params); s_fin], [vec(points); point(d.path, s_fin)]
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
