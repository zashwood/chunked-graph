using LightGraphs

import Base: sizehint!, show

mutable struct MultiGraph
	graph::Graph{Int}
	vertex_map::Dict{Label, Int}
	inverse_vertex_map::Dict{Int, Label}
	edge_map::Dict{Tuple{Int, Int}, Set{AtomicEdge}}

	function MultiGraph()
		return new(Graph(), Dict{Label, Int}(), Dict{Int, Label}(), Dict{Tuple{Int, Int}, Set{AtomicEdge}}())
	end
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

@inline function add_vertex!(mgraph::MultiGraph, lbl::Label)
	LightGraphs.add_vertex!(mgraph.graph)
	vcnt = nv(mgraph.graph)
	mgraph.vertex_map[lbl] = vcnt
	mgraph.inverse_vertex_map[vcnt] = lbl
end

@inline function rem_vertex!(mgraph::MultiGraph, lbl::Label)
	u = mgraph.vertex_map[lbl]
	for v in collect(neighbors(mgraph.graph, u))
		LightGraphs.rem_edge!(mgraph.graph, u, v)
		delete!(mgraph.edge_map, unordered(u, v))
	end
	delete!(mgraph.vertex_map, lbl)
	delete!(mgraph.inverse_vertex_map, u)
end

@inline function add_edge!(mgraph::MultiGraph, lbl_u::Label, lbl_v::Label, e::AtomicEdge)
	u = mgraph.vertex_map[lbl_u]
	v = mgraph.vertex_map[lbl_v]
	LightGraphs.add_edge!(mgraph.graph, u, v)
	uv = unordered(u, v)
	if !haskey(mgraph.edge_map, uv)
		mgraph.edge_map[uv] = Set{AtomicEdge}()
	end
	push!(mgraph.edge_map[uv], e)
end

@inline function add_edges!(mgraph::MultiGraph, lbl_u::Label, lbl_v::Label, edges::Vector{AtomicEdge})
	u = mgraph.vertex_map[lbl_u]
	v = mgraph.vertex_map[lbl_v]
	LightGraphs.add_edge!(mgraph.graph, u, v)
	uv = unordered(u, v)
	if !haskey(mgraph.edge_map, uv)
		mgraph.edge_map[uv] = Set{AtomicEdge}()
	end
	union!(mgraph.edge_map[uv], edges)
end

@inline function rem_edge!(mgraph::MultiGraph, lbl_u::Label, lbl_v::Label, e::AtomicEdge)
	u = mgraph.vertex_map[lbl_u]
	v = mgraph.vertex_map[lbl_v]
	if has_edge(mgraph.graph, u, v)
		uv = unordered(u, v)
		delete!(mgraph.edge_map[uv], e)
		if length(mgraph.edge_map[uv]) === 0
			LightGraphs.rem_edge!(mgraph.graph, u, v)
			delete!(mgraph.edge_map, uv)
		end
	end
end

# TODO: Check and simplify
function incident_edges(mgraph::MultiGraph, lbl::Label)
	u = mgraph.vertex_map[lbl]
	return chain([mgraph.edge_map[unordered(u, v)] for v in neighbors(mgraph.graph, u)]...)
end

# TODO: Check and simplify
"Returns a dictionary from edges to list of atomic edges"
function induced_edges(mgraph::MultiGraph, lbls::Vector{Label})
	vertices = Int[mgraph.vertex_map[lbl] for lbl in lbls if haskey(mgraph.vertex_map, lbl)]
	vertex_set = Set{Int}(vertices)

	ret = Dict{Tuple{Label, Label}, Set{AtomicEdge}}()
	for u in vertex_set
		for v in neighbors(mgraph.graph, u)
			if u < v && v in vertex_set
				edge = (mgraph.inverse_vertex_map[u], mgraph.inverse_vertex_map[v])
				ret[edge] = mgraph.edge_map[(u, v)]
			end
		end
	end
	return ret
end

# TODO: Check and simplify
function connected_components(mgraph::MultiGraph, lbls::Set{Label})
	lgraph = mgraph.graph
	vertices = map(x->mgraph.vertex_map[x], lbls)
	visited = Set{Label}()
	sizehint!(visited, length(vertices))
	components = Vector{Int}[]
	to_visit = Set{Label}()

	for v in vertices
		if !(v in visited)
			next_component = Label[]
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
	return Vector{Label}[map(x->mgraph.inverse_vertex_map[x], y) for y in components]
end
