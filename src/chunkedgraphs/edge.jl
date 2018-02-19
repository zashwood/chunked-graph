import Base: isless, isequal, ==, <, hash, collect, union!, setdiff!, isempty, convert

struct AtomicEdge
	u::Label
	v::Label
	affinity::Affinity
	AtomicEdge(u::Label, v::Label, affinity::Affinity) = u < v ? new(u, v, affinity) : new(v, u, affinity)
end

#An EdgeSet represents a set of AtomicEdges between leaves of u and v
struct CompositeEdgeSet
	u::Label
	v::Label
	max_affinity::Affinity
	children::Array{CompositeEdgeSet}
	nonempty::Bool
	#this should probably have size at most 4 from geometric constraints.
	#todo: we can probably keep this sorted to speed up union! and setdiff!
end

function set_nonempty(e::CompositeEdgeSet, b::Bool)
	return CompositeEdgeSet(e.u, e.v, e.max_affinity,e.children,b)
end

const EMPTY_EDGE_LIST = Vector{CompositeEdgeSet}()

function isleaf(e::CompositeEdgeSet)
	@assert (tolevel(e) == 1) == (e.children === EMPTY_EDGE_LIST)
	return tolevel(e) == 1
end

function tolevel(e::CompositeEdgeSet)
	return tolevel(head(e))
end

function Base.convert(::Type{AtomicEdge}, e::CompositeEdgeSet)
	@assert isleaf(e)
	@assert !isempty(e)
	return AtomicEdge(e.u,e.v,e.max_affinity)
end

#TODO: don't allocate so many intermediate arrays
function collect(c::CompositeEdgeSet)
	if isleaf(c)
		return [convert(AtomicEdge, c)]
	else
		return vcat(map(collect, c.children)...)::Array{AtomicEdge}
	end
end

function CompositeEdge(c::ChunkedGraph, ae::AtomicEdge)
	@assert hasvertex!(c, head(ae))
	@assert hasvertex!(c, tail(ae))
	e=CompositeEdgeSet(head(ae),tail(ae), ae.affinity, EMPTY_EDGE_LIST, true)
	while !isvalid(c, e)
		@assert hasvertex!(c,head(e))
		@assert hasvertex!(c,tail(e))
		e=CompositeEdgeSet(force_get_parent!(c,head(e)), 
						force_get_parent!(c,tail(e)),
						e.max_affinity,
						CompositeEdgeSet[e],
						true
						)
	end
	return e
end

function union!(e1::CompositeEdgeSet, e2::CompositeEdgeSet)
	@assert isleaf(e1) == isleaf(e2)
	@assert head(e1)==head(e2)
	@assert tail(e1)==tail(e2)
	if isleaf(e1)
		e1.isempty = e1.isempty || e2.isempty
	else
		for c2 in e2.children
			found=false
			for c1 in e1.children
				if (head(c1),tail(c1)) == (head(c2),tail(c2))
					union!(c1,c2)
					found=true
					break
				end
			end
			if !found
				push!(e1.children, c2)
			end
		end
	end
	return e1
end

function isempty(e::CompositeEdgeSet)
	if isleaf(e)
		return !e.nonempty
	else
		return length(e.children) == 0
	end
end

function setdiff!(e1::CompositeEdgeSet, e2::CompositeEdgeSet)
	@assert isleaf(e1) == isleaf(e2)
	@assert head(e1)==head(e2)
	@assert tail(e1)==tail(e2)

	if isleaf(e1) && e2.nonempty
		return set_nonempty(e1, false)
	else
		for c2 in e2.children
			for (i,c1) in enumerate(e1.children)
				if (head(c1),tail(c1)) == (head(c2),tail(c2))
					e1.children[i]=setdiff!(c1,c2)
					break
				end
			end
		end
		filter!(!isempty, e1.children)
	end
	return e1
end

"returns a set of validated edges equivalent to e"
function revalidate!(c::ChunkedGraph, e::CompositeEdgeSet)
	return collect(map(x->buildup!(c,x),breakdown!(c,e)))
end

function buildup!(c::ChunkedGraph, e::CompositeEdgeSet)
	while !isvalid(c, e)
		@assert hasvertex!(c, head(e))
		@assert hasvertex!(c, tail(e))
		e=CompositeEdgeSet(
					force_get_parent!(c,head(e)), 
					force_get_parent!(c,tail(e)),
					e.max_affinity,
					CompositeEdgeSet[e],
					true
					)
	end
	return e
end
function breakdown!(c::ChunkedGraph, e::CompositeEdgeSet)
	if !hasvertex!(c, head(e)) || !hasvertex!(c, tail(e)) 
		return chain(map(x->breakdown!(c,x),e.children)...)
	else
		return [e]
	end
end

head(e::CompositeEdgeSet)=e.u
tail(e::CompositeEdgeSet)=e.v
head(e::AtomicEdge)=e.u
tail(e::AtomicEdge)=e.v

function isvalid(cgraph, e)
	return 	hasvertex!(cgraph, head(e)) && 
			hasvertex!(cgraph, tail(e)) &&
			(tochunkid(head(e)) != tochunkid(tail(e)) || tolevel(head(e)) == 1) &&
			parent(tochunkid(head(e))) == parent(tochunkid(tail(e)))
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

#hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, hash(:AtomicEdge, seed)))
hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, seed))
