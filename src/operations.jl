# These definitions enable basic arithmetic operations on rational functions. Currently, they
# do not apply to approximations on domains, because testing the equality of domains is tricky.

# Promote to the widest numeric type:
Base.promote_rule(::Type{Barycentric{S}}, ::Type{Barycentric{T}}) where {T,S} = Barycentric{promote_type(S,T)}
Base.promote_rule(::Type{T}, ::Type{Barycentric{S}}) where {T<:Number,S} = Barycentric{promote_type(S,T)}

# Convert the numeric type:
function Base.convert(::Type{Barycentric{T}}, r::Barycentric{S}) where {S,T}
    return Barycentric{T}( T.(nodes(r)), T.(values(r)), T.(weights(r)), T.(r.w_times_f) )
end

# Create a constant rational function:
function Base.convert(::Type{Barycentric{T}}, x::Number) where T
    S = promote_type(T, typeof(x))
    return Barycentric{S}([zero(S)], [S(x)], [one(S)], [S(x)])
end

# Unary minus is simple:
Base.:-(r::Barycentric{T}) where T = Barycentric{T}(nodes(r), -values(r), weights(r), -r.w_times_f)

# Arithmetic with constant functions is usually also simple:
Base.:+(r::Barycentric, x::Number) = Barycentric(nodes(r), values(r) .+ x, weights(r))
Base.:+(x::Number, r::Barycentric) = r + x
Base.:-(r::Barycentric, x::Number) = Barycentric(nodes(r), values(r) .- x, weights(r))
Base.:-(x::Number, r::Barycentric) = (-r) + x
Base.:*(r::Barycentric, x::Number) = Barycentric(nodes(r), values(r) .* x, weights(r))
Base.:*(x::Number, r::Barycentric) = r * x
Base.:/(r::Barycentric, x::Number) = Barycentric(nodes(r), values(r) ./ x, weights(r))

# Other arithmetic operations require resampling:
Base.:/(x::Number, r::Barycentric) = barycentric( z -> x / r(z) )
Base.:+(r::Barycentric, s::Barycentric) = barycentric( z -> r(z) + s(z) )
Base.:-(r::Barycentric, s::Barycentric) = barycentric( z -> r(z) - s(z) )
Base.:*(r::Barycentric, s::Barycentric) = barycentric( z -> r(z) * s(z) )
Base.:/(r::Barycentric, s::Barycentric) = barycentric( z -> r(z) / s(z) )
