import Base: isequal, ==, hash

mutable struct Vertex
	label::Label
	parent::Label
	children::Vector{Label}
	#internal_edges::Array{Edge}
	#invariant: all edges of internal_edges should be maximal edges
	#invariant: children should form a connected graph via edges
end

# Comparison
==(lhs::Vertex, rhs::Vertex) = isequal(lhs.label, rhs.label)
isequal(lhs::Vertex, rhs::Vertex) = lhs == rhs

hash(v::Vertex, seed::UInt) = hash(v.label, hash(:Vertex, seed))

# Conversion / Getters
@inline function tochunk(v::Vertex)
	return tochunk(v.label)
end

@inline function tosegment(v::Vertex)
	return tosegment(v.label)
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