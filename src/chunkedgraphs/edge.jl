import Base: isless, isequal, ==, <, hash, collect, union!, setdiff!

struct AtomicEdge
	u::Label
	v::Label
	affinity::Affinity
end

type SingletonEdgeSet
	e::AtomicEdge
	nonempty::Bool
end

#An Edge represents a set of AtomicEdges between leaves of u and v
type CompositeEdgeSet
	u::Label
	v::Label
	children::Array{Union{CompositeEdgeSet,SingletonEdgeSet}} #this should probably have size at most 4 from geometric constraints.
end

const EdgeSet = Union{SingletonEdgeSet, CompositeEdgeSet}

#=
function SingletonEdgeSet()
	SingletonEdgeSet(AtomicEdge(NULL_LABEL, NULL_LABEL,1f0),true)
end
function EdgeSet()
	return CompositeEdgeSet(NULL_LABEL, NULL_LABEL, EdgeSet[])
end
=#



function collect(c::CompositeEdgeSet)
	return cat(1, map(collect, c.children)...)::Array{AtomicEdge}
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
	while !is_valid(c, e)
		e=CompositeEdgeSet(force_get_parent!(c,head(e)), 
						force_get_parent!(c,tail(e)),EdgeSet[e])
	end
	return e
end

function force_get_parent!(c::ChunkedGraph, l::Label)
	v = getvertex!(c,l)
	if v.parent == NULL_LABEL
		promote!(c, v)
	end
	return v.parent
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

	if tolevel(head(e1)) > 1
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
	else
		union!(e1.children, e2.children)
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
	while !is_valid(c, e)
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
	@assert is_valid(c,e)
	return e
end

head(e::CompositeEdgeSet)=e.u
tail(e::CompositeEdgeSet)=e.v
head(e::SingletonEdgeSet)=head(e.e)
tail(e::SingletonEdgeSet)=tail(e.e)
head(e::AtomicEdge)=e.u
tail(e::AtomicEdge)=e.v

function is_valid(chunked_graph, e)
	return 	hasvertex!(chunked_graph, head(e)) && 
			hasvertex!(chunked_graph, tail(e)) &&
			(tochunk(head(e)) != tochunk(tail(e)) || tolevel(head(e)) == 1) &&
			parent(tochunk(head(e))) == parent(tochunk(tail(e)))
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