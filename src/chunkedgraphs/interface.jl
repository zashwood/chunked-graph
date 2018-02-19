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
	c = getchunk!(cgraph, parent(tochunkid(vertex)))
	@assert haskey(c.vertices, vertex.label)
	@assert tolevel(c) < MAX_DEPTH

	@assert c.clean
	@assert vertex.parent == NULL_LABEL
	@assert length(incident_edges(c.graph, vertex.label)) == 0

	l = uniquelabel!(c)
	pv = Vertex(l, NULL_LABEL, Label[vertex.label])
	vertex.parent = pv.label

	@assert tochunkid(pv) == c.id
	add_vertex!(c.parent.graph, pv.label)
	c.parent.vertices[pv.label] = pv
	c.modified=true
	@assert hasvertex!(cgraph, vertex.parent)
	@assert tochunkid(vertex.parent) == parent(tochunkid(vertex.label))
	return pv
end

function force_get_parent!(c::ChunkedGraph, l::Label)
	v = getvertex!(c,l)
	if v.parent == NULL_LABEL
		promote!(c, v)
	end
	@assert hasvertex!(c, v.parent)
	return v.parent
end

#=
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
=#

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
	@assert tolevel(tochunkid(lbl)) == 1 "Vertex label at level $(tolevel(tochunkid(lbl))), expected 1."

	c = getchunk!(cgraph, parent(tochunkid(lbl)))
	if haskey(c.vertices, lbl)
		#TODO: warn user
		return
	end

	v = Vertex(lbl, NULL_LABEL, EMPTY_LABEL_LIST)
	add_vertex!(c, v)
end

function add_atomic_vertices!(cgraph::ChunkedGraph, lbls::Vector{Label})
	gc_enable(false)
	for (i,lbl) in enumerate(lbls)
		if i%100000==0
			println(i)
		end
		if length(cgraph.chunks) > 2*CACHESIZE
			run_eviction!(cgraph)
			gc_enable(true)
			gc()
			gc_enable(false)
		end
		add_atomic_vertex!(cgraph, lbl)
	end
	gc_enable(true)
end

function add_atomic_edge!(cgraph::ChunkedGraph, edge::AtomicEdge)
	cid = lca(tochunkid(edge.u), tochunkid(edge.v))
	if tolevel(cid) == 1
		cid = parent(cid)
	end
	c = getchunk!(cgraph, cid)
	add_edge!(c, edge)
end

function add_atomic_edges!(cgraph::ChunkedGraph, edges::Vector{AtomicEdge})
	gc_enable(false)
	for (i,edge) in enumerate(edges)
		if i%100000==0
			println(i)
		end
		if length(cgraph.chunks) > 2*CACHESIZE
			run_eviction!(cgraph)
			gc_enable(true)
			gc()
			gc_enable(false)
		end
		add_atomic_edge!(cgraph, edge)
	end
	gc_enable(true)
end

function add_atomic_edges!(cgraph::ChunkedGraph, edges::Vector{Tuple{Label,Label}})
	gc_enable(false)
	for (i,edge) in enumerate(edges)
		if i%100000==0
			println(i)
		end
		if length(cgraph.chunks) > 2*CACHESIZE
			run_eviction!(cgraph)
			gc_enable(true)
			gc()
			gc_enable(false)
		end
		add_atomic_edge!(cgraph, AtomicEdge(edge[1],edge[2]))
	end
	gc_enable(true)
end

function delete_atomic_edge!(cgraph::ChunkedGraph, edge::AtomicEdge)
	cid = lca(tochunkid(edge.u), tochunkid(edge.v))
	if tolevel(cid) == 1
		cid = parent(cid)
	end
	c = getchunk!(cgraph, cid)
	delete_edge!(c, edge)
end
