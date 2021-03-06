function leaves!(cgraph::ChunkedGraph, vertex::Vertex, stop_lvl::Integer = 1, bbox::Union{Cuboid, Void} = nothing)
	@assert tolevel(tochunkid(vertex)) >= stop_lvl

	if tolevel(tochunkid(vertex)) == stop_lvl
		return [vertex.label]
	end

	gc_enable(false)
	vertices = [vertex]
	lvl = tolevel(tochunkid(vertex))
	while lvl > stop_lvl + 1
		vertices = map(lbl->getvertex!(cgraph, lbl), mapreduce(v->v.children, vcat, Label[], vertices))
		if bbox !== nothing
			filter!(v->overlaps(tocuboid(tochunkid(v)), bbox::Cuboid), vertices)
		end
		lvl -= 1
	end
	ret = mapreduce(v->v.children, vcat, Label[], vertices)

	gc_enable(true)
	return ret
end

function promote!(cgraph::ChunkedGraph, vertex::Vertex)
	c = getchunk!(cgraph, tochunkid(vertex))
	@assert tolevel(c) < MAX_DEPTH

	@assert c.clean
	@assert vertex.parent == NULL_LABEL
	@assert length(incident_edges(c.graph, vertex.label)) == 0

	l = uniquelabel!(c.parent)
	pv = Vertex(l, NULL_LABEL, Label[vertex.label])
	vertex.parent = pv.label

	@assert tochunkid(pv) == c.parent.id
	add_vertex!(c.parent.graph, pv.label)
	c.parent.vertices[pv.label] = pv
	c.parent.modified = true

	return pv
end

function promote_to_lca!(cgraph::ChunkedGraph, vertex1::Vertex, vertex2::Vertex)
	if tochunkid(vertex1) == tochunkid(vertex2)
		return (vertex1, vertex2)
	else
		if vertex1.parent == NULL_LABEL
			promote!(cgraph, vertex1)
		end
		if vertex2.parent == NULL_LABEL
			promote!(cgraph, vertex2)
		end
		parent1 = getvertex!(cgraph, vertex1.parent)
		parent2 = getvertex!(cgraph, vertex2.parent)
		return promote_to_lca!(cgraph, parent1, parent2)
	end
end

function root!(cgraph::ChunkedGraph, vertex::Vertex)
	if isroot(vertex)
		return vertex
	else
		return root!(cgraph, getvertex!(cgraph, vertex.parent))
	end
end

n_processed = 0
function update!(cgraph::ChunkedGraph)
	global n_processed
	n_processed = 0
	gc_enable(false)
	update!(getchunk!(cgraph, TOP_ID))
	gc_enable(true)
end

function add_atomic_vertex!(cgraph::ChunkedGraph, lbl::Label)
	@assert tolevel(tochunkid(lbl)) == 1 "Vertex label at level $(tolevel(tochunkid(label))), expected 1."

	c = getchunk!(cgraph, tochunkid(lbl))
	if haskey(c.vertices, lbl)
		#TODO: warn user
		return
	end

	v = Vertex(lbl, NULL_LABEL, EMPTY_LABEL_LIST)
	push!(c.added_vertices, v)
	touch!(c)
end

function add_atomic_vertices!(cgraph::ChunkedGraph, lbls::Vector{Label})
	gc_enable(false)
	for lbl in lbls
		add_atomic_vertex!(cgraph, lbl)
	end
	gc_enable(true)
end

function add_atomic_edge!(cgraph::ChunkedGraph, edge::AtomicEdge)
	c = getchunk!(cgraph, lca(tochunkid(edge.u), tochunkid(edge.v)))
	push!(c.added_edges, edge)
	touch!(c)
end

function add_atomic_edges!(cgraph::ChunkedGraph, edges::Vector{AtomicEdge})
	gc_enable(false)
	for edge in edges
		add_atomic_edge!(cgraph, edge)
	end
	gc_enable(true)
end

function delete_atomic_edge!(cgraph::ChunkedGraph, edge::AtomicEdge)
	c = getchunk!(cgraph, lca(tochunkid(edge.u), tochunkid(edge.v)))
	push!(c.deleted_edges, edge)
	touch!(c)
end
