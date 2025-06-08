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

### DiscretizedPath

struct DiscretizedPath{T,F}
    path::T              # original Curve or Path
    points::Matrix       # first column = active points, others are refinements
    params::Matrix{F}    # first column = active points, others are refinements
    next::Vector{Int}    # next[i] = index of the next point in the path
end
"""
    DiscretizedPath(path, s::AbstractVector; kwargs...)
    DiscretizedPath(path, n::Integer=0; kwargs...)

Discretize a path, keeping the option of future making local refinements.

# Arguments
- `path`: a ComplexCurve or ComplexPath
- `s`: a vector of parameter values
- `n`: number of points to discretize the path

# Keyword arguments
- `refinement`: number of refinements to make between consecutive points
- `maxpoints`: maximum number of points ever allowed

See also [`collect`](@ref), [`add_node!`](@ref).
"""
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

function DiscretizedPath(p::DiscretizedPath, refinement::Int)
    s, _ = collect(p, :nodes)
    maxpoints = max(length(s), size(p.points, 1))
    return DiscretizedPath(p.path, s; refinement, maxpoints)
end
"""
    add_node!(d::DiscretizedPath, idx)

Add a new node to the discretization, and return the indexes of all affected points. The indexes
are valid on the `points` and `params` fields.

# Arguments
- `d`: a DiscretizedPath object
- `idx`: a 2-element tuple, vector, or `CartesianIndex` into the `params` field. This identifies
the point to be promoted to a node.

# Returns
- A 2-element vector of `CartesianIndices` into the `params` and `points` fields.

See also [`DiscretizedPath`](@ref), [`collect`](@ref).
"""
function add_node!(d::DiscretizedPath, idx)
    if length(idx) != 2
        throw(ArgumentError("Second argument must be a 2-element tuple or CartesianIndex"))
    end
    (idx[2] == 1) && return
    n = length(d.next) + 1
    if n == size(d.points, 1)
        throw(error("Cannot add more points to this discretization"))
    end
    if idx[2] > size(d.params, 2)
        throw(error("Indicated new node is not in the discretization"))
    end

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
end

"""
    collect(d::DiscretizedPath, which=:nodes)

Collect the points and parameters of a discretized path.

# Arguments
- `d`: a DiscretizedPath object
- `which`: return the nodes if :nodes, test points if :test, or all if :all

# Returns
- Tuple of two vectors: parameter values and points on the path

See also [`add_node!`](@ref), [`DiscretizedPath`](@ref).
"""
function Base.collect(d::DiscretizedPath, which=:nodes)
    n = length(d.next)
    if which == :nodes
        columns = [1]
    elseif which == :test
        columns = 2:size(d.params, 2)
    elseif which == :all
        columns = 1:size(d.params, 2)
    else
        throw(ArgumentError("Second argument must be :nodes, :test, or :all"))
    end
    params = similar(d.params, length(columns), n)
    points = similar(d.points, length(columns), n)
    k = idx = 1
    while idx != 0
        params[:, k] .= d.params[idx, columns]
        points[:, k] .= d.points[idx, columns]
        idx = d.next[idx]
        k += 1
    end
    if which in (:nodes, :all)
        s_fin = length(d.path)
        return [vec(params); s_fin], [vec(points); point(d.path, s_fin)]
    else
        return vec(params), vec(points)
    end
end

# Find the distance to the nearest neighbor at every point.
function spacing(d::DiscretizedPath{T,F}) where {T,F}
    m = length(d.next)
    n = size(d.params, 2)
    Δ = fill(F(Inf), m, n)
    k = idx = 1
    while true
        δ = @. abs(d.points[idx, 2:n] - d.points[idx, 1:n-1])
        @. Δ[idx, 2:n-1] = min(δ[1:n-2], δ[2:n-1])
        Δ[idx, 1] = min(Δ[idx, 1], δ[1])
        idx1 = d.next[idx]
        if idx1 == 0
            Δ[idx, n] = δ[n-1]
            break
        end
        δ₊ = abs(d.points[idx1, 1] - d.points[idx, n])
        Δ[idx, n] = min(δ[n-1], δ₊)
        Δ[idx1, 1] = δ₊
        idx = idx1
        k += 1
    end
    return Δ
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
