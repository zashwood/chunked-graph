import Base: collect, union!, setdiff!

const SingletonEdgeSet = Union{AtomicEdge, Void}

#An EdgeSet represents a set of AtomicEdges between leaves of u and v
struct CompositeEdgeSet
	u::Label
	v::Label
	children::Vector{Union{SingletonEdgeSet, CompositeEdgeSet}}
	#this should probably have size at most 4 from geometric constraints.
	#todo: we can probably keep this sorted to speed up union! and setdiff!
end

const EdgeSet = Union{SingletonEdgeSet, CompositeEdgeSet}

collect(e::SingletonEdgeSet) = typeof(e) === Void ? AtomicEdge[] : AtomicEdge[e]
function collect(c::CompositeEdgeSet)
	#TODO: don't allocate so many intermediate arrays
	return vcat(map(collect, c.children)...)::Vector{AtomicEdge}
end

function CompositeEdge(c::ChunkedGraph, e::AtomicEdge)
	@assert hasvertex!(c, e.u)
	@assert hasvertex!(c, e.v)
	while !isvalid(c, e)
		e = CompositeEdgeSet(force_get_parent!(c, e.u), force_get_parent!(c, e.v), EdgeSet[e])
	end
	return e
end

union!(e1::SingletonEdgeSet, e2::Void) = e1
union!(e1::SingletonEdgeSet, e2::AtomicEdge) = e2
function union!(e1::CompositeEdgeSet, e2::CompositeEdgeSet)
	@assert e1.u === e2.u
	@assert e1.v === e2.v
	@assert tolevel(e1.u) > 1
	for c2 in e2.children
		found = false
		for (i, c1) in enumerate(e1.children)
			if (c1.u, c1.v) === (c2.u, c2.v)
				e1.children[i] = union!(c1, c2)
				found = true
				break
			end
		end
		if !found
			push!(e1.children, c2)
		end
	end
	return e1
end

setdiff!(e1::SingletonEdgeSet, e2::Void) = e1
setdiff!(e1::SingletonEdgeSet, e2::AtomicEdge) = nothing
function setdiff!(e1::CompositeEdgeSet, e2::CompositeEdgeSet)
	@assert e1.u === e2.u
	@assert e1.v === e2.v

	for c2 in e2.children
		for (i, c1) in enumerate(e1.children)
			if (c1.u, c1.v) === (c2.u, c2.v)
				e1.children[i] = setdiff!(c1, c2)
				break
			end
		end
	end
	filter!(!isempty, e1.children)
	return e1
end

isempty(e::CompositeEdgeSet) = length(e.children) === 0
isempty(e::SingletonEdgeSet) = typeof(e) === Void

"returns a set of validated edges equivalent to e"
function revalidate!(c::ChunkedGraph, e::CompositeEdgeSet)
	return collect(map(x->buildup!(c, x), breakdown!(c, e)))
end
function revalidate!(c::ChunkedGraph, e::AtomicEdge)
	@assert isvalid(c, e)
	return e
end

function buildup!(c::ChunkedGraph, e::EdgeSet)
	@assert hasvertex!(c, e.u)
	@assert hasvertex!(c, e.v)
	while !isvalid(c, e)
		e = CompositeEdgeSet(force_get_parent!(c, e.u),
		                     force_get_parent!(c, e.v), EdgeSet[e])
	end
	return e
end

function breakdown!(c::ChunkedGraph, e::CompositeEdgeSet)
	if !hasvertex!(c, e.u) || !hasvertex!(c, e.v)
		return chain(map(x->breakdown!(c, x), e.children)...)
	else
		return [e]
	end
end

function breakdown!(c::ChunkedGraph, e::AtomicEdge)
	@assert hasvertex!(c, e.u) && hasvertex!(c, e.v)
	return [e]
end

function isvalid(cgraph::ChunkedGraph, e::EdgeSet)
	return hasvertex!(cgraph, e.u) && hasvertex!(cgraph, e.v) &&
	       (tochunkid(e.u) !== tochunkid(e.v) || tolevel(e.u) == 1) &&
	       parent(tochunkid(e.u)) === parent(tochunkid(e.v))
end
