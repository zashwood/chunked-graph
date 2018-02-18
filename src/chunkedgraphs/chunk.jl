using IterTools

mutable struct Chunk
	cgraph::ChunkedGraph{Chunk}
	id::ChunkID
	graph::MultiGraph{Label, AtomicEdge}
	vertices::Dict{Label, Vertex}
	parent::Union{Chunk, Void}
	children::Vector{Chunk}
	added_vertices::Set{Vertex}
	deleted_vertices::Set{Vertex}
	added_edges::Set{AtomicEdge}
	deleted_edges::Set{AtomicEdge}
	clean::Bool
	modified::Bool
	max_label::SegmentID

	#todo: lazy loading for vertices

	function Chunk(cgraph::ChunkedGraph, chunkid::ChunkID, vertices::Dict{Label,Vertex}, graph::MultiGraph, max_label::Label)
		c = new(cgraph, chunkid, graph, vertices, nothing, Chunk[], Set{Vertex}(), Set{Vertex}(), Set{AtomicEdge}(), Set{AtomicEdge}(), true, false, max_label)

		if !isroot(chunkid)
			par = parent!(c)
			push!(par.children, c)
			c.parent = par
		end
		return c
	end
end

@inline function add_vertex!(c::Chunk, v::Vertex)
	@assert !(v in c.deleted_vertices)
	@assert parent(tochunkid(v.label)) === c.id
	push!(c.added_vertices, v)
end

@inline function delete_vertex!(c::Chunk, v::Vertex)
	@assert !(v in c.added_vertices)
	@assert parent(tochunkid(v.label)) === c.id
	push!(c.deleted_vertices, v)
end

@inline function add_edge!(c::Chunk, e::AtomicEdge)
	@assert !(e in c.deleted_edges)
	@assert begin
		tmp = lca(tochunkid(e.u), tochunkid(e.v))
		if tolevel(tmp) == 1
			tmp = parent(tmp)
		end
		tmp == c.id
	end
	push!(c.added_edges, e)
end

@inline function delete_edge!(c::Chunk, e::AtomicEdge)
	@assert !(e in c.added_edges)
	@assert begin
		tmp = lca(tochunkid(e.u), tochunkid(e.v))
		if tolevel(tmp) == 1
			tmp = parent(tmp)
		end
		tmp == c.id
	end
	push!(c.deleted_edges, e)
end

function Chunk(cgraph::ChunkedGraph, chunkid::ChunkID)
	m = MultiGraph{Label, AtomicEdge, tolevel(chunkid) == 1 ? SingletonEdgeSet : CompositeEdgeSet}()
	return Chunk(cgraph, chunkid, Dict{Label, Vertex}(), m, Label(0))
end

@inline function tochunkid(c::Chunk)
	return c.id
end

@inline function tolevel(c::Chunk)
	return tolevel(c.id)
end

@inline function topos(c::Chunk)
	return topos(c.id)
end

@inline function isroot(c::Chunk)
	return isroot(c.id)
end

@inline function parent!(c::Chunk)
	return getchunk!(c.cgraph, parent(c.id))
end

function touch!(c::Void)
end
function touch!(c::Chunk)
	if c.clean
		c.clean = false
		c.modified = true
		touch!(c.parent)
	end
end

function uniquelabel!(c::Chunk)
	#TODO: Should erase clean flag?
	c.max_label += 1
	return tolabel(c.id, c.max_label)
end

#TODO: use a queue to maintain order of operations
function update!(c::Chunk)
	#TODO: tag updates with time since updates don't commute
	if c.clean
		return
	end

	# Update children first
	for child in c.children
		update!(child)
	end

	# Print debug messages
	# println("updating $(tolevel(c.id)), $(map(Int,topos(c.id))) V: +$(length(c.added_vertices))/-$(length(c.deleted_vertices)), E: +$(length(c.added_edges))/-$(length(c.deleted_edges))")
	global n_processed
	n_processed += 1
	# println("$n_processed/$(length(c.cgraph.chunks))")

	#vertices which need updates
	dirty_vertices = Set{Vertex}()
	redo_edge_sets = Set{EdgeSet}()

	# FIXME: We should upsize with the difference of added minus deleted
	upsize!(c.graph, length(c.added_vertices), length(c.added_edges))
	upsize!(c.vertices, length(c.added_vertices))

	# Insert added vertices
	# mark them as dirty_vertices as well
	for v in c.added_vertices
		@assert parent(tochunkid(v)) == c.id
		@assert v.parent == NULL_LABEL
		@assert !haskey(c.vertices, v.label)
		add_vertex!(c.graph, v.label)
		c.vertices[v.label] = v
		push!(dirty_vertices,v)
	end

	# Delete vertices and mark the vertices connected to 
	# the one we are deleting as dirty
	for v in c.deleted_vertices
		# for child in v.children
		# 	@assert get_vertex(c.cgraph, child).parent === NULL_LABEL
		# end

		for edge_set in incident_edges(c.graph, v.label)
			push!(redo_edge_sets, edge_set)
		end

		# this should delete all edges incident with v as well
		rem_vertex!(c.graph, v.label)
		delete!(c.vertices, v.label)
	end

	for e in redo_edge_sets
		for ee in revalidate!(c.cgraph, e)
			u = c.vertices[ee.u]
			v = c.vertices[ee.v]
			add_edges!(c.graph, u.label, v.label, ee)
			push!(dirty_vertices, u, v)
		end
	end

	for e in c.added_edges
		e = CompositeEdge(c.cgraph, e)
		u, v = c.vertices[e.u], c.vertices[e.v]
		@assert tochunkid(u.label) !== tochunkid(v.label) || (tolevel(u.label) == 1 && tolevel(v.label) == 1)
		@assert parent(tochunkid(u.label)) === parent(tochunkid(v.label)) === tochunkid(c)
		@assert haskey(c.vertices, u.label)
		@assert haskey(c.vertices, v.label)

		@assert haskey(c.graph.vertex_map, u.label)
		@assert haskey(c.graph.vertex_map, v.label)
		add_edges!(c.graph, u.label, v.label, e)
		push!(dirty_vertices, u, v)
	end

	for e in c.deleted_edges
		e = CompositeEdge(c.cgraph, e)
		u, v = c.vertices[e.u], c.vertices[e.v]
		@assert isvalid(c.cgraph, e)
		@assert parent(tochunkid(e.u)) === parent(tochunkid(e.v)) === tochunkid(c)
		@assert tochunkid(u.label) !== tochunkid(v.label) || (tolevel(u.label) == 1 && tolevel(v.label) == 1)
		@assert parent(tochunkid(u.label)) === parent(tochunkid(v.label)) === tochunkid(c)
		rem_edges!(c.graph, u.label, v.label, e)
		push!(dirty_vertices, u, v)
	end

	if !isroot(c)
		for v in chain(dirty_vertices, c.deleted_vertices)
			if v.parent !== NULL_LABEL
				@assert tochunkid(v.parent) === tochunkid(c)
				@assert parent(tochunkid(v.parent)) === tochunkid(c.parent)
				delete_vertex!(c.parent, c.parent.vertices[v.parent])
				v.parent = NULL_LABEL
			end
		end
		c.parent.clean = false
	end

	cc = connected_components(c.graph, map(x->x.label, dirty_vertices))
	for component in cc
		if length(component) > 1
			l = uniquelabel!(c)
			new_vertex = Vertex(l, NULL_LABEL, component)
			for child_label in component
				c.vertices[child_label].parent = new_vertex.label
			end
			if !isroot(c)
				push!(c.parent.added_vertices, new_vertex)
			end
		else
			c.vertices[component[1]].parent = NULL_LABEL
		end
	end

	for v in values(c.vertices)
		@assert v.parent === NULL_LABEL || tochunkid(v.parent) === c.id
	end

	empty!(c.added_edges)
	empty!(c.added_vertices)
	empty!(c.deleted_edges)
	empty!(c.deleted_vertices)
	c.clean = true
end
