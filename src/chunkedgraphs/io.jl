import Base: show

show(io::IO, cgraph::ChunkedGraph) = print("ChunkedGraph with $(length(cgraph.chunks)) in memory")

function write(io::IO, e::CompositeEdgeSet)
	write(io, e.u)
	write(io,e.v)
	write(io,length(e.children))
	for c in e.children
		write(io,c)
	end
end

function read(io::IO, e::CompositeEdgeSet)
	u=read(io,Label)
	v=read(io,Label)
	n=read(io,Int64)
	T=(tolevel(u)>=2 ? CompositeEdgeSet : SingletonEdgeSet)
	children = T[read(io, T) for i in 1:n]
	return CompositeEdgeSet(u, v, children)
end

function write(io::IO, e::SingletonEdgeSet)
	write(io,e.u)
	write(io,e.v)
	write(io,e.nonempty)
end

function read(io::IO, e::SingletonEdgeSet)
	return SingletonEdgeSet(read(io,Label),read(io,Label),read(io,Bool))
end

function loadchunk(cgraph::ChunkedGraph, chunkid::ChunkID)
	vertex_map = Dict{Label, Vertex}()
	mgraph=MultiGraph{Label,AtomicEdge,tolevel(chunkid)==2 ? SingletonEdgeSet : CompositeEdgeSet}()
	max_label = Label(0)

	prefix = stringify(chunkid)
	path = expanduser(joinpath(cgraph.path, String("$(prefix).chunk")))
	if !isfile(path)
		return Chunk(cgraph, chunkid, vertex_map, mgraph, max_label)
	end

	
	f = open(path, "r")

	# Check File Version
	version = read(f, UInt64)
	@assert version === UInt64(1)

	# Read Chunk Info
	(max_label, v_cnt, e_cnt) = read(f, UInt64, 3)

	# Allocate Graph
	sizehint!(mgraph, floor(UInt32, 1.5 * v_cnt), floor(UInt32, 1.5 * e_cnt))
	sizehint!(vertex_map, floor(UInt32, 1.5 * v_cnt))
	
	# Read Vertices
	for i in range(1, v_cnt)
		(label, parent, child_cnt) = read(f, UInt64, 3)
		children = read(f, UInt64, child_cnt)
		add_vertex!(mgraph, label)
		vertex_map[label] = Vertex(label, parent, children)
	end
	T=tolevel(chunkid) == 2 ? SingletonEdgeSet : CompositeEdgeSet
	# Read EdgeMap
	for i in range(1, e_cnt)
		edge_set = read(f, T)
		add_edges!(mgraph, e.u, e.v, edge_set)
	end

	close(f)
	return Chunk(cgraph, chunkid, vertex_map, mgraph, max_label)
end

function getchunk!(cgraph::ChunkedGraph, chunkid::ChunkID)
	global eviction_mode
	if !haskey(cgraph.chunks, chunkid)
		tmp = loadchunk(cgraph, chunkid)
		if tmp == nothing
			tmp = Chunk(cgraph, chunkid)
		end
		cgraph.chunks[chunkid] = tmp::Chunk
	end

	c = cgraph.chunks[chunkid]
	if !eviction_mode
		if tolevel(c) == 1
			cgraph.lastused[c.id] = time()
		end
		
		eviction_mode = true
		while length(cgraph.chunks) > CACHESIZE
			evict!(cgraph.chunks[DataStructures.peek(cgraph.lastused)[1]])
		end
		eviction_mode = false
	else
		if tolevel(c) == 1 && !haskey(cgraph.lastused, c.id)
			cgraph.lastused[c.id] = time()
		end
	end
	return c
end

function getvertex!(cgraph::ChunkedGraph, l::Label)
	return getchunk!(cgraph, parent(tochunkid(l))).vertices[l]
end

function hasvertex!(cgraph::ChunkedGraph, l::Label)
	return haskey(getchunk!(cgraph, parent(tochunkid(l))).vertices,l)
end

function save!(cgraph::ChunkedGraph, force::Bool = false)
	for c in collect(values(cgraph.chunks))
		if c.modified || force
			save!(c)
		end
	end
end

function save!(c::Chunk)
	@assert c.clean

	prefix = stringify(c.id)
	path = expanduser(joinpath(c.cgraph.path, "$(prefix).chunk"))
	print("Saving to $(path)...")

	buf = IOBuffer()
	write(buf, UInt64(1)) # Version
	write(buf, UInt64(c.max_label)) # Max SegID
	write(buf, UInt64(length(c.vertices))) # Vertex Count
	write(buf, UInt64(length(LightGraphs.edges(c.graph.g)))) # Edge Count

	for vertex in values(c.vertices)
		write(buf, UInt64(vertex.label)) # Vertex Label
		write(buf, UInt64(vertex.parent)) # Vertex Parent
		write(buf, UInt64(length(vertex.children))) # Vertex Children Count
		write(buf, convert(Vector{UInt64}, vertex.children))
	end
		
	for edgeset in values(c.graph.edge_map)
		write(buf, edgeset)
	end

	f = open(path, "w")
	write(f, buf.data)
	close(f)
	close(buf)
	println("done.")
	
	c.modified = false
end

function evict!(c::Chunk)
	@assert !isroot(c)
	@assert c.id != SECOND_ID
	update!(c)
	while !isempty(c.children)
		evict!(c.children[end])
	end

	if c.modified
		save!(c)
	end

	filter!(x->x != c, c.parent.children)
	delete!(c.cgraph.chunks, c.id)
	dequeue!(c.cgraph.lastused, c.id)
end
