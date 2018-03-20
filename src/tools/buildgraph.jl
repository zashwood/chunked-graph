push!(LOAD_PATH, dirname(@__FILE__))
include("../chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs
using CloudVolume

import PyCall
import JSON

rel(p::String) = joinpath(dirname(@__FILE__), p)
settings = JSON.parsefile(rel("graph.conf"))

# Ensure graph directory exists and is empty
if !isdir(settings["graphpath"])
	mkdir(settings["graphpath"])
end
if !isempty(readdir(settings["graphpath"]))
	error("'$(settings["graphpath"])' is not empty!")
end

# Init Cloud Storage and Chunked Graph
edgetasks = StorageWrapper(settings["edgetasks"])
graph = ChunkedGraph(settings["graphpath"], "gs://neuroglancer/removeme/wow")

# Build Graph
ranges = ([settings["offset"][i] : settings["step"][i]: settings["offset"][i] + settings["size"][i] for i in 1:3]...)
total = length(ranges[1]) * length(ranges[2]) * length(ranges[3])
@time begin
	i = 0
	for z in ranges[3], y in ranges[2], x in ranges[1]
		i += 1
		prefix = "$(x)-$(x+settings["step"][1])_$(y)-$(y+settings["step"][2])_$(z)-$(z+settings["step"][3])"
		atomicedges = edgetasks.val[:get_file]("$(prefix)_atomicedges.bin")
		#rg2cg = edgetasks.val[:get_file]("$(prefix)_rg2cg.bin")

		atomicedges = Vector{AtomicEdge}(reinterpret(AtomicEdge, Vector{UInt8}(atomicedges)))
		println("$i/$total | $prefix: Adding $(length(atomicedges)) edges")
		#rg2cg = Dict{UInt64, Label}(reinterpret(Pair{UInt64, Label}, Vector{UInt8}(rg2cg)))
		
		add_atomic_vertices!(graph, unique(vcat(map(x->x.u, atomicedges), map(x->x.v, atomicedges))))
		add_atomic_edges!(graph, atomicedges)
	end
	update!(graph)
	save!(graph)
end


#DEBUGGING...

# evict! = ChunkedGraphs.evict!
# upsize! = ChunkedGraphs.upsize!
# add_vertex! = ChunkedGraphs.add_vertex!
# rem_vertex! = ChunkedGraphs.rem_vertex!
# add_edges! = ChunkedGraphs.add_edges!
# revalidate! = ChunkedGraphs.revalidate!
# incident_edges = ChunkedGraphs.incident_edges
# parent = ChunkedGraphs.parent
# head = ChunkedGraphs.head
# tail = ChunkedGraphs.tail
# CACHESIZE = 25
# const eviction_count = DataStructures.DefaultDict{ChunkID,Int}(0)
# cgraph = graph

# gc_enable(true)
# cgraph.eviction_mode = true
# while length(cgraph.chunks) > CACHESIZE
# 	convict, priority = DataStructures.peek(cgraph.lastused)
# 	children = cgraph.chunks[convict].children
# 	if length(children) > 0
# 		#spare him, kill the children first
# 		#he'll be put back on death row when his children die
# 		DataStructures.dequeue!(cgraph.lastused, convict)
# 		continue
# 	else
# 		evict!(cgraph.chunks[convict])
# 		if (length(cgraph.chunks) == 34)
# 			break
# 		end
# 	end
# end

# convict, priority = DataStructures.peek(cgraph.lastused)
# children = cgraph.chunks[convict].children;

# c = cgraph.chunks[convict];
# t = (eviction_count[c.id] += 1)

# Vertex = ChunkedGraphs.Vertex
# CompositeEdgeSet = ChunkedGraphs.CompositeEdgeSet
# dirty_vertices = Set{Vertex}()
# redo_edge_sets = Set{CompositeEdgeSet}()

# # FIXME: We should upsize with the difference of added minus deleted
# upsize!(c.graph, length(c.added_vertices), length(c.added_edges))
# upsize!(c.vertices, length(c.added_vertices))

# # Insert added vertices
# # mark them as dirty_vertices as well
# for v in c.added_vertices
# 	@assert parent(tochunkid(v)) == c.id
# 	@assert v.parent == NULL_LABEL
# 	@assert !haskey(c.vertices, v.label)
# 	add_vertex!(c.graph, v.label)
# 	c.vertices[v.label] = v
# 	push!(dirty_vertices,v)
# end

# # Delete vertices and mark the vertices connected to 
# # the one we are deleting as dirty
# for v in c.deleted_vertices
# 	# for child in v.children
# 	# 	@assert get_vertex(c.cgraph, child).parent === NULL_LABEL
# 	# end

# 	for edge_set in incident_edges(c.graph, v.label)
# 		push!(redo_edge_sets, edge_set)
# 	end

# 	# this should delete all edges incident with v as well
# 	rem_vertex!(c.graph, v.label)
# 	delete!(c.vertices, v.label)
# end

# redo_edge = 0
# for e in redo_edge_sets
# 	if e.u == 0x0601000000000001 && e.v == 0x06010100000000a7
# 		redo_edge = e
# 		break
# 	end
# 	for ee in revalidate!(c.cgraph, e)
# 		u = c.vertices[head(ee)]
# 		v = c.vertices[tail(ee)]
# 		add_edges!(c.graph,u.label,v.label,ee)
# 		push!(dirty_vertices, u, v)
# 	end
# end

# isvalid = ChunkedGraphs.isvalid
# broken = ChunkedGraphs.breakdown!(c.cgraph, redo_edge)
# bad_edge = first(broken)
# map(x -> ChunkedGraphs.buildup2!(c.cgraph, x), broken)

# # while !isempty(c.children)
# # 	ChunkedGraphs.evict!(c.children[end])
# # end

# # if c.modified
# # 	ChunkedGraphs.save!(c)
# # end

# # filter!(x->x != c, c.parent.children)
# # delete!(c.cgraph.chunks, c.id)
# # priority = c.cgraph.lastused[c.id]
# # DataStructures.dequeue!(c.cgraph.lastused, c.id)
# # if length(c.parent.children) == 0 && !haskey(c.cgraph.lastused, c.parent.id)
# # 	c.cgraph.lastused[c.parent.id] = priority+1
# # end
# # unlock(c.flock)
# # close(c.flock)