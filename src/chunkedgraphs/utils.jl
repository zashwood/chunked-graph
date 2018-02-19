@inline function unordered(uv::Tuple{T, T}) where T <: Integer
	@inbounds return minmax(uv[1], uv[2])
end

@inline function unordered(u::T, v::T) where T <: Integer
	return minmax(u, v)
end

@inline function upsize!(a::AbstractArray, n::Integer)
	sizehint!(a, length(a) + n)
end

@inline function upsize!(a::Associative, n::Integer)
	sizehint!(a, length(a) + n)
end
