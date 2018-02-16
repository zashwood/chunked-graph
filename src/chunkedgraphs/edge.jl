import Base: isless, isequal, ==, <, hash, collect, union!, setdiff!

struct AtomicEdge
	u::Label
	v::Label
	affinity::Affinity
	AtomicEdge(u::Label, v::Label, affinity::Affinity) = u < v ? new(u, v, affinity) : new(v, u, affinity)
end

mutable struct SingletonEdgeSet
	e::AtomicEdge
	nonempty::Bool
end

#An EdgeSet represents a set of AtomicEdges between leaves of u and v
struct CompositeEdgeSet
	u::Label
	v::Label
	children::Union{Array{CompositeEdgeSet},Array{SingletonEdgeSet}} 
	#this should probably have size at most 4 from geometric constraints.
	#todo: we can probably keep this sorted to speed up union! and setdiff!
end

const EdgeSet = Union{SingletonEdgeSet, CompositeEdgeSet}

#TODO: don't allocate so many intermediate arrays
function collect(c::CompositeEdgeSet)
	return vcat(map(collect, c.children)...)::Array{AtomicEdge}
end

function collect(c::SingletonEdgeSet)
	if c.nonempty
		return AtomicEdge[c.e]
	else
		return AtomicEdge[]
	end
end

function CompositeEdge(c::ChunkedGraph, e::AtomicEdge)
	@assert hasvertex!(c, head(e))
	@assert hasvertex!(c, tail(e))
	e=SingletonEdgeSet(e,true)
	while !isvalid(c, e)
		e=CompositeEdgeSet(force_get_parent!(c,head(e)), 
						force_get_parent!(c,tail(e)),EdgeSet[e])
	end
	return e
end

function union!(e1::SingletonEdgeSet, e2::SingletonEdgeSet)
	@assert head(e1)==head(e2)
	@assert tail(e1)==tail(e2)
	if !isempty(e2)
		e1.nonempty=true
	end

	return e1
end

function setdiff!(e1::SingletonEdgeSet, e2::SingletonEdgeSet)
	@assert head(e1)==head(e2)
	@assert tail(e1)==tail(e2)
	if !isempty(e2)
		e1.nonempty=false	
	end
	return e1
end

function union!(e1::CompositeEdgeSet, e2::CompositeEdgeSet)
	@assert head(e1)==head(e2)
	@assert tail(e1)==tail(e2)
	@assert tolevel(head(e1)) > 1
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
	return e1
end

function isempty(e::CompositeEdgeSet)
	return length(e.children) == 0
end
function isempty(e::SingletonEdgeSet)
	return !e.nonempty
end

function setdiff!(e1::CompositeEdgeSet, e2::CompositeEdgeSet)
	@assert head(e1)==head(e2)
	@assert tail(e1)==tail(e2)

	for c2 in e2.children
		for c1 in e1.children
			if (head(c1),tail(c1)) == (head(c2),tail(c2))
				setdiff!(c1,c2)
				break
			end
		end
	end
	filter!(!isempty, e1.children)
	return e1
end

"returns a set of validated edges equivalent to e"
function revalidate!(c::ChunkedGraph, e::CompositeEdgeSet)
	return collect(map(x->buildup!(c,x),breakdown!(c,e)))
end

function buildup!(c::ChunkedGraph, e::Union{CompositeEdgeSet,SingletonEdgeSet})
	@assert hasvertex!(c, head(e))
	@assert hasvertex!(c, tail(e))
	while !isvalid(c, e)
		e=CompositeEdgeSet(force_get_parent!(c,head(e)), 
						force_get_parent!(c,tail(e)),EdgeSet[e])
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
function breakdown!(c::ChunkedGraph, e::SingletonEdgeSet)
	@assert hasvertex!(c, head(e)) && hasvertex!(c, tail(e))
	return [e]
end
function revalidate!(c::ChunkedGraph, e::SingletonEdgeSet)
	@assert isvalid(c,e)
	return e
end

head(e::CompositeEdgeSet)=e.u
tail(e::CompositeEdgeSet)=e.v
head(e::SingletonEdgeSet)=head(e.e)
tail(e::SingletonEdgeSet)=tail(e.e)
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

hash(e::AtomicEdge, seed::UInt) = hash(e.u, hash(e.v, hash(:AtomicEdge, seed)))