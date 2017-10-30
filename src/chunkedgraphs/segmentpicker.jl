function getsupervoxelat(cgraph::ChunkedGraph, start_label::Label, voxel_pos::Tuple{Int, Int, Int})
	if tochunk(start_label) == 0
		start_label = tolabel(world_to_chunk(voxel_pos[1], voxel_pos[2], voxel_pos[3]), tosegment(start_label))
		# FIXME: If we deal with very coarse supervoxel meshes near chunk boundaries, the picked position might lie within
		#        the neighboring chunk which doesn't contain the supervoxel
	end

	const root_vertex = root!(cgraph, getvertex!(cgraph, start_label))

	const padding = CartesianIndex(16,16,2)
	const voxel_res = Vector{Int}(cgraph.cloudvolume.val[:scale]["resolution"])
	const cutout_range = voxel_pos[1]-padding[1] : voxel_pos[1]+padding[1], voxel_pos[2]-padding[2] : voxel_pos[2]+padding[2], voxel_pos[3]-padding[3] : voxel_pos[3]+padding[3]

	const segmentation_ids = cgraph.cloudvolume[cutout_range...]::Array{UInt32, 3} # <-- slow, even when cached
	const world_positions = collect(CartesianRange(cutout_range))

	const sortkey = sortperm(collect(Base.Iterators.flatten(
		map(i->(voxel_res[1] * (i[1] - padding[1] - 1))^2 + (voxel_res[2] * (i[2] - padding[2] - 1))^2 + (voxel_res[3] * (i[3] - padding[3] - 1))^2,
				CartesianRange(size(segmentation_ids))
		)
	)))

	const already_checked = Set{Int}()
	for i in sortkey
		if segmentation_ids[i] == 0 || segmentation_ids[i] in already_checked
			continue
		end

		segid = segmentation_ids[i]
		worldpos = world_positions[i]
		chunkid = world_to_chunk(worldpos[1], worldpos[2], worldpos[3])
		labelid = tolabel(chunkid, segid)

		if root!(cgraph, getvertex!(cgraph, labelid)) == root_vertex
			return labelid
		end

		push!(already_checked, segid)
	end

	return nothing
end