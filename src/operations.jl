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
