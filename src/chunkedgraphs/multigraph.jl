using LightGraphs
#using SmallSets

import Base: sizehint!, show

#C is a container type for elements of type E
mutable struct MultiGraph{L,E,C}
	graph::Graph{Int}
	vertex_map::Dict{L, Int}
	inverse_vertex_map::Dict{Int, L}
	edge_map::Dict{Tuple{Int, Int}, C}
end

function MultiGraph{L,E,C}() where {L,E,C}
	return MultiGraph{L,E,C}(Graph(), Dict{L, Int}(), Dict{Int, L}(), Dict{Tuple{Int, Int}, C}())
end

show(io::IO, mgraph::MultiGraph) = print("MultiGraph with $(nv(mgraph.graph)) and $(ne(mgraph.graph)) edges")

@inline function sizehint!(mgraph::MultiGraph, v_cnt::Integer, e_cnt::Integer)
	sizehint!(mgraph.vertex_map, v_cnt)
	sizehint!(mgraph.inverse_vertex_map, v_cnt)
	sizehint!(mgraph.graph.fadjlist, v_cnt)
	sizehint!(mgraph.edge_map, e_cnt)	
end

@inline function upsize!(mgraph::MultiGraph, v_cnt::Integer, e_cnt::Integer)
	upsize!(mgraph.vertex_map, v_cnt)
	upsize!(mgraph.inverse_vertex_map, v_cnt)
	upsize!(mgraph.graph.fadjlist, v_cnt)
	upsize!(mgraph.edge_map, e_cnt)
end

@inline function add_vertex!(mgraph::MultiGraph{L,E}, lbl::L) where {L,E}
	LightGraphs.add_vertex!(mgraph.graph)
	vcnt = nv(mgraph.graph)
	mgraph.vertex_map[lbl] = vcnt
	mgraph.inverse_vertex_map[vcnt] = lbl
end

@inline function rem_vertex!(mgraph::MultiGraph{L,E}, lbl::L) where {L,E}
	u = mgraph.vertex_map[lbl]
	for v in collect(neighbors(mgraph.graph, u))
		LightGraphs.rem_edge!(mgraph.graph, u, v)
		delete!(mgraph.edge_map, unordered(u, v))
	end
	delete!(mgraph.vertex_map, lbl)
	delete!(mgraph.inverse_vertex_map, u)
end

function add_edge!(mgraph::MultiGraph{L,E,C}, lbl_u::L, lbl_v::L, e::E) where {L,E,C}
	u = mgraph.vertex_map[lbl_u]
	v = mgraph.vertex_map[lbl_v]
	LightGraphs.add_edge!(mgraph.graph, u, v)
	uv = unordered(u, v)
	if !haskey(mgraph.edge_map, uv)
		mgraph.edge_map[uv] = C()
	end
	push!(mgraph.edge_map[uv], e)
end

function add_edges!(mgraph::MultiGraph{L,E,C}, lbl_u::L, lbl_v::L, edges::Union{Set{E},Vector{E},C}) where {L,E,C}
	u = mgraph.vertex_map[lbl_u]
	v = mgraph.vertex_map[lbl_v]
	LightGraphs.add_edge!(mgraph.graph, u, v)
	uv = unordered(u, v)
	if !haskey(mgraph.edge_map, uv)
		mgraph.edge_map[uv] = edges
	else
		union!(mgraph.edge_map[uv], edges)
	end
end

function rem_edge!(mgraph::MultiGraph{L,E}, lbl_u::L, lbl_v::L, e::E) where {L,E}
	u = mgraph.vertex_map[lbl_u]
	v = mgraph.vertex_map[lbl_v]
	if has_edge(mgraph.graph, u, v)
		uv = unordered(u, v)
		delete!(mgraph.edge_map[uv], e)
		if isempty(mgraph.edge_map[uv])
			LightGraphs.rem_edge!(mgraph.graph, u, v)
			delete!(mgraph.edge_map, uv)
		end
	end
end

function rem_edges!(mgraph::MultiGraph{L,E,C}, lbl_u::L, lbl_v::L, edges::Union{Set{E},Vector{E},C}) where {L,E,C}
	u = mgraph.vertex_map[lbl_u]
	v = mgraph.vertex_map[lbl_v]
	if has_edge(mgraph.graph, u, v)
		uv = unordered(u, v)
		setdiff!(mgraph.edge_map[uv], edges)
		if isempty(mgraph.edge_map[uv])
			LightGraphs.rem_edge!(mgraph.graph, u, v)
			delete!(mgraph.edge_map, uv)
		end
	end
end

# TODO: Check and simplify
function incident_edges(mgraph::MultiGraph{L,E}, lbl::L) where {L,E}
	u = mgraph.vertex_map[lbl]
	return [mgraph.edge_map[unordered(u, v)] for v in neighbors(mgraph.graph, u)]
end

# TODO: Check and simplify
"Returns a list of EdgeSets among vertices of lbls"
function induced_edges(mgraph::MultiGraph{L,E,C}, lbls::Vector{L}) where {L,E,C}
	vertices = Int[mgraph.vertex_map[lbl] for lbl in lbls if haskey(mgraph.vertex_map, lbl)]
	vertex_set = Set{Int}(vertices)

	ret = C[]
	for u in vertex_set
		for v in neighbors(mgraph.graph, u)
			if u < v && v in vertex_set
				edge = (mgraph.inverse_vertex_map[u], mgraph.inverse_vertex_map[v])
				push!(ret, mgraph.edge_map[(u, v)])
			end
		end
	end
	return ret
end

"Returns the atomic edges induced by lbls"
function induced_atomic_edges(mgraph::MultiGraph{L,E}, lbls::Vector{L}) where {L,E}
	return cat(1, (collect(i) for i in induced_edges(mgraph,lbls))...)
end


# TODO: Check and simplify
function connected_components(mgraph::MultiGraph{L,E}, lbls::Set{L}) where {L,E}
	lgraph = mgraph.graph
	vertices = map(x->mgraph.vertex_map[x], lbls)
	visited = Set{L}()
	sizehint!(visited, length(vertices))
	components = Vector{Int}[]
	to_visit = Set{L}()

	for v in vertices
		if !(v in visited)
			next_component = L[]
			empty!(to_visit)
			push!(to_visit, v)

			while length(to_visit) > 0
				x = pop!(to_visit)
				push!(next_component, x)
				push!(visited, x)
				for n in neighbors(lgraph, x)
					if !(n in visited)
						push!(to_visit, n)
					end
				end
			end
			push!(components, next_component)
		end
	end
	#@assert length(vertices) == sum(map(length,components))
	return Vector{L}[map(x->mgraph.inverse_vertex_map[x], y) for y in components]
end