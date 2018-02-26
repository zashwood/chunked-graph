push!(LOAD_PATH, dirname(@__FILE__))
include("../chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs
using CloudVolume
using CodecZstd
using Memoize
using MySQL
using OffsetArrays

import PyCall

mutable struct RegionGraphWrapper
	storage::StorageWrapper
	chunksize::Tuple{Unsigned, Unsigned, Unsigned}
	offset::Tuple{Integer, Integer, Integer}
	maxlevel::Unsigned
end

const RG_CHUNKSIZE = (512, 512, 512)
const CG_CHUNKSIZE = (512, 512, 64)
const mysql_conn = MySQL.connect("35.231.65.155", "root", "**********"; db="relabeling")

@inline function drop_last(x::UnitRange{T}) where T <: Integer
	return x.start : x.stop - 1
end

@inline function stringify(x::UnitRange{T}) where T <: Integer
	return "$(x.start)-$(x.stop)"
end

@inline function stringify(x::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T <: Integer
	return join(map(stringify, x), '_')
end

@inline function toslice(s::AbstractString)
	m = match(r"(\d+)-(\d+)_(\d+)-(\d+)_(\d+)-(\d+)", s)
	bounds = map(x -> parse(Int, x), m.captures)
	return (bounds[1] : bounds[2], bounds[3] : bounds[4], bounds[5] : bounds[6])
end

function load_labels(watershed::OffsetArray{T, 3}) where T <: Integer
	voxel_range = indices(watershed)[1:3]
	chunk_range = map(x -> div.(x, CG_CHUNKSIZE), (minimum.(voxel_range), maximum.(voxel_range)))
	res = MySQL.query(mysql_conn, """
		SELECT c.x, c.y, c.z, m.regiongraph_sv, m.chunkedgraph_sv
		FROM chunks c JOIN mappings m ON c.id = m.chunk_id
		WHERE c.x BETWEEN $(chunk_range[1][1]) AND $(chunk_range[2][1] + 1) AND
		      c.y BETWEEN $(chunk_range[1][2]) AND $(chunk_range[2][2] + 1) AND
			  c.z BETWEEN $(chunk_range[1][3]) AND $(chunk_range[2][3] + 1);
	""", DataFrame)

	chunk_labels = Dict{ChunkedGraphs.ChunkID, Dict{T, Label}}()
	for row in eachrow(res)
		chunkid = tochunkid(1, row[:x], row[:y], row[:z])



	end
	return chunk_labels::Dict{ChunkedGraphs.ChunkID, Dict{T, Label}}
end

function save_labels(chunk_labels::Dict{ChunkedGraphs.ChunkID, Dict{T, Label}})
	
end

function chunked_labelling(watershed::OffsetArray{T, 3}) where T <: Integer
	# load existing labels of neighboring + center chunk

	relabeled = OffsetArray(Label, indices(watershed)...)
	max_id = 1 # TODO: Currently fixed to UInt32
	rg_to_cg = Dict{Tuple{ChunkedGraphs.ChunkID, T}, Label}()
	cg_to_rg = Dict{Label, T}()
	for k in indices(watershed, 3), j in indices(watershed, 2), i in indices(watershed, 1)
		chunkid = ChunkedGraphs.world_to_chunk(i, j, k)
		if watershed[i, j, k] == 0
			# Don't relabel boundary
			relabeled[i, j, k] = 0
		elseif haskey(rg_to_cg, watershed[i, j, k])
			relabeled[i, j, k] = rg_to_cg[(chunkid, watershed[i, j, k])]
		else
			rg_to_cg[(chunkid, watershed[i, j, k])] = relabeled[i, j, k] = tolabel(chunkid, UInt32(max_id))
			max_id += 1 # Depends on chunk!
		end
		cg_to_rg[relabeled[i, j, k]] = watershed[i, j, k]
	end

	# save new labels of neighboring + center chunk
	return relabeled, cg_to_rg, rg_to_cg
end

function compute_regiongraph(watershed::AbstractArray{S, 3}, relabeled::AbstractArray{Label, 3}, agglomeration::AbstractArray{T, 3},
		regiongraph_edges::Set{AtomicEdge}) where {S <: Integer, T <: Integer}
	edges = Set{AtomicEdge}()
	@inbounds @simd for k in drop_last(indices(relabeled, 3))
		for j in indices(relabeled, 2), i in indices(relabeled, 1)
			if relabeled[i, j, k] !== relabeled[i, j, k + 1] && agglomeration[i, j, k] === agglomeration[i, j, k + 1]
				if watershed[i, j, k] === watershed[i, j, k + 1]
					push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i, j, k + 1], INF_CAPACITY))
				else
					edge = getkey(regiongraph_edges.dict, AtomicEdge(watershed[i, j, k], watershed[i, j, k + 1]), nothing)
					if edge isa AtomicEdge
						push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i, j, k + 1], edge.affinity))
					else
						println("Ruh roh for $(AtomicEdge(watershed[i, j, k], watershed[i, j, k + 1]))")
					end
				end
			end
		end
	end
	@inbounds @simd for k in indices(relabeled, 3)
		for j in drop_last(indices(relabeled, 2)), i in indices(relabeled, 1)
			if relabeled[i, j, k] !== relabeled[i, j + 1, k] && agglomeration[i, j, k] === agglomeration[i, j + 1, k]
				if watershed[i, j, k] === watershed[i, j + 1, k]
					push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i, j + 1, k], INF_CAPACITY))
				else
					edge = getkey(regiongraph_edges.dict, AtomicEdge(watershed[i, j, k], watershed[i, j + 1, k]), nothing)
					if edge isa AtomicEdge
						push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i, j + 1, k], edge.affinity))
					else
						println("Ruh roh for $(AtomicEdge(watershed[i, j, k], watershed[i, j + 1, k]))")
					end
				end
			end
		end
	end
	@inbounds @simd for k in indices(relabeled, 3)
		for j in indices(relabeled, 2), i in drop_last(indices(relabeled, 1))
			if relabeled[i, j, k] !== relabeled[i + 1, j, k] && agglomeration[i, j, k] === agglomeration[i + 1, j, k]
				if watershed[i, j, k] === watershed[i + 1, j, k]
					push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i + 1, j, k], INF_CAPACITY))
				else
					edge = getkey(regiongraph_edges.dict, AtomicEdge(watershed[i, j, k], watershed[i + 1, j, k]), nothing)
					if edge isa AtomicEdge
						push!(edges, AtomicEdge(relabeled[i, j, k], relabeled[i + 1, j, k], edge.affinity))
					else
						println("Ruh roh for $(AtomicEdge(watershed[i, j, k], watershed[i + 1, j, k]))")
					end
				end
			end
		end
	end
	return collect(edges)
end

# struct RanStruct
# 	segA1::UInt64
# 	segB1::UInt64
# 	sum_aff1::Float32
# 	sum_area1::UInt64
# 	segA2::UInt64
# 	segB2::UInt64
# 	sum_aff2::Float32
# 	sum_area2::UInt64
# end
#
# The first 4 values are currently the same as the last 4 values...
# Affinity between segA1 and segB1 is sum_aff1 / sum_area1
# Values were written without padding, which is why we can't read them in Julia as Array

@memoize Dict function get_chunk_affinities(storage::StorageWrapper, chunk_path::AbstractString)
	edges = Vector{AtomicEdge}()

	f = storage.val[:get_file](chunk_path)
	if f isa Void
		warn("$chunk_path doesn't exist")
		return edges
	end

	buf = ZstdDecompressorStream(IOBuffer(f))
	sizehint!(edges, 500000)
	while !eof(buf)
		segA, segB = read(buf, UInt64, 2); aff = read(buf, Float32); area = read(buf, UInt64); skip(buf, 28);
		push!(edges, AtomicEdge(segA, segB, aff / area))
	end

	return edges
end

function get_chunkhierarchy_affinities(regiongraph::RegionGraphWrapper, slices::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T
	chunk_range = map(x -> div.(x .- regiongraph.offset, regiongraph.chunksize), (minimum.(slices), maximum.(slices)))

	edges = Vector{AtomicEdge}()
	sizehint!(edges, 1000000)
	for l in 0:regiongraph.maxlevel
		for x in chunk_range[1][1]:chunk_range[2][1], y in chunk_range[1][2]:chunk_range[2][2], z in chunk_range[1][3]:chunk_range[2][3]
			chunk_path = "complete_edges_$(l)_$(x)_$(y)_$(z).data.zst"
			println(chunk_path)
			append!(edges, get_chunk_affinities(regiongraph.storage, chunk_path))
		end
		chunk_range = (div.(first(chunk_range), 2), div.(last(chunk_range), 2))
	end

	return edges
end

function edge_task(watershed::CloudVolumeWrapper, agglomeration::CloudVolumeWrapper, regiongraph::RegionGraphWrapper,
		output_storage::StorageWrapper, slices::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T<:Integer

	# no need to fetch _inner_ edges from the neighboring chunks - the boundary crossing edges are stored in higher levels
	regiongraph_edges = get_chunkhierarchy_affinities(regiongraph, slices)

	watershed_cutout = OffsetArray(watershed[slices...], slices...)
	agglomeration_cutout = OffsetArray(agglomeration[slices...], slices...)

	relabeled_cutout, cg_to_rg, rg_to_cg = chunked_labelling(watershed_cutout)
	edges = compute_regiongraph(watershed_cutout, relabeled_cutout, agglomeration_cutout, regiongraph_edges)
	vertices = unique(relabeled_cutout)

	edge_buf = IOBuffer()
	vertex_buf = IOBuffer()
	write(edge_buf, edges)
	write(vertex_buf, vertices)
	output_storage.val[:put_file](file_path="$(stringify(slices))_edges.bin", content = PyCall.pybytes(edge_buf.data))
	output_storage.val[:put_file](file_path="$(stringify(slices))_vertices.bin", content = PyCall.pybytes(vertex_buf.data))
	output_storage.val[:wait]()
end

function edge_task_bundle(watershed_path::AbstractString, agglomeration_path::AbstractString,
		regiongraph_path::AbstractString, output_path::AbstractString,
		slices::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where S <: Unsigned where T <: Integer
	slices = map(drop_last, slices) # julia indexing is inclusive

	watershed = CloudVolumeWrapper(watershed_path; bounded=false)
	agglomeration = CloudVolumeWrapper(agglomeration_path; bounded=false)

	regiongraph = RegionGraphWrapper(
		StorageWrapper(regiongraph_path),
		RG_CHUNKSIZE,
		(offset(watershed)...),
		Int(ceil(log2(maximum( div.(size(watershed)[1:3], RG_CHUNKSIZE) )))) # max octree level
	)
	output_storage = StorageWrapper(output_path)

	for x_start in first(slices[1]) : CG_CHUNKSIZE[1] : last(slices[1])
		for y_start in first(slices[2]) : CG_CHUNKSIZE[2] : last(slices[2])
			for z_start in first(slices[3]) : CG_CHUNKSIZE[3] : last(slices[3])
				subslice = (x_start : x_start + CG_CHUNKSIZE[1],
				            y_start : y_start + CG_CHUNKSIZE[2],
				            z_start : z_start + CG_CHUNKSIZE[3])
				edge_task(watershed, agglomeration, regiongraph, output_storage, subslice)
			end
		end
	end

end

watershed_path, agglomeration_path, output_path, slice_str = ARGS

watershed_path = "gs://neuroglancer/ranl/pinky40_watershed_test"
agglomeration_path = "gs://neuroglancer/ranl/pinky40_agglomeration_test_20"
regiongraph_path = "gs://ranl/pinky40_agglomeration_test/region_graph"
output_path = "gs://nkem/pinky40_agglomeration_test/region_graph"


edge_task(watershed_path, agglomeration_path, output_path, toslice(slice_str))
