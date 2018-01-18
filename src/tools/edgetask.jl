push!(LOAD_PATH, dirname(@__FILE__))
include("../chunkedgraphs/ChunkedGraphs.jl")
using ChunkedGraphs
using CloudVolume
using OffsetArrays

import PyCall

const pyslice = PyCall.pybuiltin(:slice)

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

function chunked_labelling(raw::OffsetArray{T, 3}) where T <: Integer
	ret = OffsetArray(Label, indices(raw)...)
	for k in indices(raw, 3), j in indices(raw, 2), i in indices(raw, 1)
		ret[i, j, k] = tolabel(ChunkedGraphs.world_to_chunk(i, j, k), raw[i, j, k])
	end
	return ret
end

function compute_regiongraph(raw::AbstractArray{T, 3}, machine_labels::AbstractArray{S, 3}) where {T <: Integer, S <: Integer}
	edges = Set{AtomicEdge}()
	@inbounds @simd for k in drop_last(indices(raw, 3))
		for j in indices(raw, 2), i in indices(raw, 1)
			if raw[i, j, k] !== raw[i, j, k + 1] && machine_labels[i, j, k] === machine_labels[i, j, k + 1]
				push!(edges, AtomicEdge(raw[i, j, k], raw[i, j, k + 1]))
			end
		end
	end
	@inbounds @simd for k in indices(raw, 3)
		for j in drop_last(indices(raw, 2)), i in indices(raw, 1)
			if raw[i, j, k] !== raw[i, j + 1, k] && machine_labels[i, j, k] === machine_labels[i, j + 1, k]
				push!(edges, AtomicEdge(raw[i, j, k], raw[i, j + 1, k]))
			end
		end
	end
	@inbounds @simd for k in indices(raw, 3)
		for j in indices(raw, 2), i in drop_last(indices(raw, 1))
			if raw[i, j, k] !== raw[i + 1, j, k] && machine_labels[i, j, k] === machine_labels[i + 1, j, k]
				push!(edges, AtomicEdge(raw[i + 1, j, k], raw[i, j, k]))
			end
		end
	end
	return collect(edges)
end

function edge_task(watershed_path::AbstractString, segmentation_path::AbstractString, output_path::AbstractString,
		slices::Tuple{UnitRange{T}, UnitRange{T}, UnitRange{T}}) where T<:Integer
	slice_x, slice_y, slice_z = slices
	watershed = CloudVolumeWrapper(watershed_path; bounded=false)
	watershed_cutout = OffsetArray(watershed[slice_x, slice_y, slice_z], slice_x, slice_y, slice_z)

	segmentation = CloudVolumeWrapper(segmentation_path; bounded=false)
	segmentation_cutout = OffsetArray(segmentation[slice_x, slice_y, slice_z], slice_x, slice_y, slice_z)

	relabelled = chunked_labelling(watershed_cutout)
	edges = compute_regiongraph(relabelled, segmentation_cutout)
	vertices = unique(relabelled)
	
	output_storage = CloudVolume.StorageWrapper(output_path)
	edge_buf = IOBuffer()
	vertex_buf = IOBuffer()
	write(edge_buf, edges)
	write(vertex_buf, vertices)
	output_storage.val[:put_file](file_path="$(stringify(slices))_edges.bin", content = PyCall.pybytes(edge_buf.data))
	output_storage.val[:put_file](file_path="$(stringify(slices))_vertices.bin", content = PyCall.pybytes(vertex_buf.data))
	output_storage.val[:wait]()
end

watershed_path, segmentation_path, output_path, slice_str = ARGS
slices = map(drop_last, toslice(slice_str)) # julia indexing is inclusive

edge_task(watershed_path, segmentation_path, output_path, slices)
