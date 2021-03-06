import Base: isequal, ==, hash

mutable struct Vertex
	label::Label
	parent::Label
	children::Vector{Label}
end

# Comparison
==(lhs::Vertex, rhs::Vertex) = isequal(lhs.label, rhs.label)
isequal(lhs::Vertex, rhs::Vertex) = lhs == rhs

hash(v::Vertex, seed::UInt) = hash(v.label, hash(:Vertex, seed))

# Conversion / Getters
@inline function tochunkid(v::Vertex)
	return tochunkid(v.label)
end

@inline function tosegid(v::Vertex)
	return tosegid(v.label)
end

@inline function tolevel(v::Vertex)
	return tolevel(v.label)
end

@inline function topos(v::Vertex)
	return topos(v.label)
end

# Misc
@inline function isleaf(v::Vertex)
	return tolevel(v) == 1
end

@inline function isroot(v::Vertex)
	return v.parent == NULL_LABEL
end