import Base: isless, isequal, ==, <, hash

struct AtomicEdge
	us::Array{Label}
	vs::Array{Label}
	aff::Affinity
	function Edge(u,v)
		u,v = min(u,v), max(u,v)
		@assert tolevel(u) == 0
		@assert tolevel(v) == 0
		return new(Label[u],Label[v])
	end
end

AtomicEdge(u, v) = AtomicEdge(Label(u), Label(v), Affinity(1))
AtomicEdge(u, v, affinity) = AtomicEdge(Label(u), Label(v), Affinity(affinity))

# Comparison
(==)(lhs::AtomicEdge, rhs::AtomicEdge) = isequal(head(lhs), head(rhs)) && isequal(tail(lhs), tail(rhs))
isequal(lhs::AtomicEdge, rhs::AtomicEdge) = lhs == rhs

#TODO: Make this definition independent of bit layout!
#We mainly want to compare on chunkid first
(<)(lhs::AtomicEdge, rhs::AtomicEdge) = isequal(head(lhs), head(rhs)) ? tail(lhs) < tail(rhs) : head(lhs) < head(rhs)
isless(lhs::AtomicEdge, rhs::AtomicEdge) = lhs < rhs

hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, hash(:AtomicEdge, seed)))

function head(e::AtomicEdge)
	return e.us[end]
end
function tail(e::AtomicEdge)
	return e.vs[end]
end

function is_valid(chunked_graph, e)
	return 	hasvertex!(chunked_graph, head(e)) && 
			hasvertex!(chunked_graph, tail(e)) &&
			(chunk_id(head(e)) != chunk_id(tail(e)) || tolevel(head(e)) == 0) &&
			parent(chunk_id(head(e))) == parent(chunk_id(tail(e)))
end

function revalidate!(chunked_graph, e)
	while !hasvertex!(chunked_graph, head(e))
		@assert length(e.us) >= 1
		e.us = e.us[1:end-1]
	end
	while !hasvertex(chunked_graph, tail(e))
		@assert length(e.vs) >= 1
		e.vs = e.vs[1:end-1]
	end
	while tolevel(head(e)) < tolevel(tail(e))
		push!(e.us, getvertex!(chunked_graph, e.us[end]).parent)
	end
	while tolevel(tail(e)) < tolevel(head(e))
		push!(e.vs, getvertex!(chunked_graph, e.vs[end]).parent)
	end
	while parent(tochunk(head(e))) != parent(tochunk(tail(e)))
		push!(e.us, getvertex!(chunked_graph, e.us[end]).parent)
		push!(e.vs, getvertex!(chunked_graph, e.vs[end]).parent)
	end
end