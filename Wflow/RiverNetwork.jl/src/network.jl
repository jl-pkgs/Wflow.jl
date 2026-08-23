# Maps `EdgeConnectivity` fields to Cartesian neighbors.
const DIRS = (:ind_y_down, :ind_x_down, :ind_x_up, :ind_y_up)
const CARTESIAN_NEIGHBORS = (
    CartesianIndex(0, -1),
    CartesianIndex(-1, 0),
    CartesianIndex(1, 0),
    CartesianIndex(0, 1),
)

"""
Struct for storing 2D staggered-grid edge connectivity in the x and y directions.
Edges without neighbors use the extra index `n + 1`.
"""
@kwdef struct EdgeConnectivity
    n::Int
    ind_x_up::Vector{Int} = zeros(Int, n)
    ind_x_down::Vector{Int} = zeros(Int, n)
    ind_y_up::Vector{Int} = zeros(Int, n)
    ind_y_down::Vector{Int} = zeros(Int, n)
end

"Struct for storing network information for the land domain."
@kwdef struct NetworkLand
    modelsize::Tuple{Int, Int} = (0, 0)
    local_drain_direction::Vector{UInt8} = UInt8[]
    allocation_area_indices::Vector{Vector{Int64}} = Vector{Int}[]
    graph::SimpleDiGraph{Int} = DiGraph(0)
    streamorder::Vector{Int} = Int[]
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    land_indices::Vector{Int} = 1:length(indices)
    order::Vector{Int} = Int[]
    order_of_subdomains::Vector{Vector{Int}} = Vector{Int}[]
    order_subdomain::Vector{Vector{Int}} = Vector{Int}[]
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    river_indices::Vector{Int} = Int[]
    river_inds_excl_reservoir::Vector{Int} = Int[]
    edge_indices::EdgeConnectivity = EdgeConnectivity(; n = 1)
    subdomain_indices::Vector{Vector{Int}} = Vector{Int}[]
    upstream_nodes::Vector{Vector{Int}} = Vector{Int}[]
end

"Struct for storing network information for the river domain."
@kwdef struct NetworkRiver
    local_drain_direction::Vector{UInt8} = UInt8[]
    allocation_area_indices::Vector{Vector{Int64}} = Vector{Int}[]
    edges_at_node::EdgesAtNode = EdgesAtNode()
    graph::SimpleDiGraph{Int} = DiGraph(0)
    streamorder::Vector{Int} = Int[]
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    land_indices::Vector{Int} = Int[]
    pit_indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    nodes_at_edge::NodesAtEdge = NodesAtEdge()
    order::Vector{Int} = Int[]
    order_of_subdomains::Vector{Vector{Int}} = Vector{Int}[]
    order_subdomain::Vector{Vector{Int}} = Vector{Int}[]
    reservoir_indices::Vector{Int} = Int[]
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    subdomain_indices::Vector{Vector{Int}} = Vector{Int}[]
    upstream_nodes::Vector{Vector{Int}} = Vector{Int}[]
end

"Struct for storing network information for reservoirs."
@kwdef struct NetworkReservoir
    indices_coverage::Vector{Vector{CartesianIndex{2}}} = Vector{CartesianIndex{2}}[]
    indices_outlet::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    reverse_indices::Matrix{Int} = zeros(Int, 0, 0)
    land_indices::Vector{Int} = Int[]
    river_indices::Vector{Int} = Int[]
end

"Struct for storing network information for groundwater drains."
@kwdef struct NetworkDrain
    indices::Vector{CartesianIndex{2}} = CartesianIndex{2}[]
    reverse_indices::Matrix{Int64} = zeros(Int, 0, 0)
    land_indices::Vector{Int} = Int[]
end

"Derive active 1D↔2D indices from a gridded domain and its `nodata` value."
function active_indices(
    domain::AbstractMatrix,
    nodata,
)::Tuple{Vector{CartesianIndex{2}}, Matrix{Int}}
    indices = filter(index -> !isequal(domain[index], nodata), CartesianIndices(domain))
    reverse_indices = zeros(Int, size(domain))
    for (node_idx, cartesian_idx) in enumerate(indices)
        reverse_indices[cartesian_idx] = node_idx
    end
    return indices, reverse_indices
end

"Set pit values in a local drainage direction vector from a gridded pit mask."
function set_pit_ldd(
    pits_2d::AbstractMatrix{Bool},
    ldd::Vector{UInt8},
    indices::Vector{CartesianIndex{2}};
    pit::Integer = LDD_PIT,
)::Vector{UInt8}
    for node_idx in eachindex(indices)
        pits_2d[indices[node_idx]] && (ldd[node_idx] = UInt8(pit))
    end
    return ldd
end

"Convert gridded drainage directions to a directed acyclic graph."
function flowgraph(
    ldd::AbstractVector,
    indices::AbstractVector,
    directions::AbstractVector = PCR_DIR,
)
    graph = DiGraph(length(indices))
    for (from_node, from_index) in enumerate(indices)
        ldd_value = ldd[from_node]
        ldd_value == LDD_PIT && continue
        to_index = from_index + directions[ldd_value]
        to_node = searchsortedfirst(indices, to_index)
        if to_node > length(indices) || indices[to_node] != to_index
            @warn "Invalid drainage direction value at node `$from_node` (LDD=`$ldd_value`), assign pit value at node"
            ldd[from_node] = LDD_PIT
            continue
        end
        add_edge!(graph, from_node, to_node)
    end
    is_cyclic(graph) && error("""One or more cycles detected in flow graph.
        The provided local drainage direction map may be unsound.
        Verify that each active flow cell flows towards a pit.
        """)
    return graph, ldd
end

"Return a drainage graph and 1D local drainage directions from gridded input."
function get_drainage_network(
    ldd_2d::AbstractMatrix,
    indices::Vector{CartesianIndex{2}};
    pits_2d::Union{Nothing, AbstractMatrix{Bool}} = nothing,
)
    ldd = UInt8.(ldd_2d[indices])
    isnothing(pits_2d) || set_pit_ldd(pits_2d, ldd, indices)
    return flowgraph(ldd, indices)
end

"Initialize a land network from gridded subcatchment and drainage-direction data."
function NetworkLand(
    subcatch_2d::AbstractMatrix,
    ldd_2d::AbstractMatrix;
    active_ids = nothing,
    nodata = missing,
    pits_2d::Union{Nothing, AbstractMatrix{Bool}} = nothing,
)
    if !isnothing(active_ids)
        ids = Set(active_ids)
        subcatch_2d = map(value -> value in ids ? value : nodata, subcatch_2d)
    end
    indices, reverse_indices = active_indices(subcatch_2d, nodata)
    graph, local_drain_direction = get_drainage_network(ldd_2d, indices; pits_2d)
    order = topological_sort_by_dfs(graph)
    return NetworkLand(;
        modelsize = size(subcatch_2d),
        indices,
        reverse_indices,
        local_drain_direction,
        graph,
        order,
        streamorder = stream_order(graph, order),
    )
end

"Initialize a river network from a river mask, drainage directions and land network."
function NetworkRiver(
    river_location::AbstractArray,
    ldd_2d::AbstractMatrix,
    network::NetworkLand;
    pits_2d::Union{Nothing, AbstractMatrix{Bool}} = nothing,
)
    river_values =
        ndims(river_location) == 2 ? river_location[network.indices] : river_location
    length(river_values) == length(network.indices) ||
        throw(DimensionMismatch("river mask must match the land network"))
    land_indices = findall(!iszero, river_values)
    indices = network.indices[land_indices]
    reverse_indices = zeros(Int, network.modelsize)
    for (river_idx, cartesian_idx) in enumerate(indices)
        reverse_indices[cartesian_idx] = river_idx
    end
    graph, local_drain_direction = get_drainage_network(ldd_2d, indices; pits_2d)
    order = topological_sort_by_dfs(graph)
    return NetworkRiver(;
        indices,
        reverse_indices,
        local_drain_direction,
        graph,
        order,
        streamorder = network.streamorder[land_indices],
        land_indices,
    )
end

"Set kinematic-wave subdomains using a minimum stream order."
function network_subdomains(network::Union{NetworkLand, NetworkRiver}, min_streamorder::Int)
    pit_indices = findall(isequal(UInt8(LDD_PIT)), network.local_drain_direction)
    order, indices, topological_order = kinwave_set_subdomains(
        network.graph,
        network.order,
        pit_indices,
        network.streamorder,
        min_streamorder,
    )
    @reset network.order_of_subdomains = order
    @reset network.order_subdomain = topological_order
    @reset network.subdomain_indices = indices
    return network
end

"Initialize staggered-grid edge connectivity for a land network."
function EdgeConnectivity(network::NetworkLand)
    (; modelsize, indices, reverse_indices) = network
    edge_indices = EdgeConnectivity(; n = length(indices))
    n_rows, n_columns = modelsize
    for (node_idx, cartesian_idx) in enumerate(indices)
        for (direction_idx, neighbor) in enumerate(CARTESIAN_NEIGHBORS)
            neighbor_idx = cartesian_idx + neighbor
            direction = DIRS[direction_idx]
            valid =
                1 <= neighbor_idx[1] <= n_rows &&
                1 <= neighbor_idx[2] <= n_columns &&
                !iszero(reverse_indices[neighbor_idx])
            getfield(edge_indices, direction)[node_idx] =
                valid ? reverse_indices[neighbor_idx] : length(indices) + 1
        end
    end
    return edge_indices
end

"Filter upstream graph neighbors using a logical exclusion vector."
function filter_upstream_nodes(
    graph::SimpleDiGraph{Int},
    excluded::Vector{Bool},
)::Vector{Vector{Int}}
    upstream_nodes = Vector{Int}[]
    for node_idx in topological_sort_by_dfs(graph)
        push!(
            upstream_nodes,
            filter(index -> !excluded[index], inneighbors(graph, node_idx)),
        )
    end
    return upstream_nodes
end

"Initialize a reservoir network from outlet IDs and gridded reservoir coverage."
function NetworkReservoir(
    outlet_ids::AbstractVector,
    coverage_2d::AbstractMatrix,
    network::NetworkRiver,
)
    coverage_indices = Vector{CartesianIndex{2}}[]
    reverse_indices = zeros(Int, size(coverage_2d))
    reservoir_map = zeros(Int, length(network.indices))
    outlet_indices = CartesianIndex{2}[]
    for (river_idx, cartesian_idx) in enumerate(network.indices)
        reservoir_id = outlet_ids[river_idx]
        reservoir_id > 0 || continue
        push!(outlet_indices, cartesian_idx)
        reservoir_idx = length(outlet_indices)
        reservoir_map[river_idx] = reservoir_idx
        reverse_indices[cartesian_idx] = reservoir_idx
        push!(coverage_indices, findall(isequal(reservoir_id), coverage_2d))
    end
    river_indices = findall(!iszero, reservoir_map)
    reservoir = NetworkReservoir(;
        indices_outlet = outlet_indices,
        indices_coverage = coverage_indices,
        reverse_indices,
        river_indices,
        land_indices = network.land_indices[river_indices],
    )
    return reservoir, reservoir_map
end

"Initialize a drain network from a land-domain drain mask."
function NetworkDrain(
    drain::AbstractVector,
    indices::Vector{CartesianIndex{2}},
    surface_flow_width::AbstractVector{<:Real},
    modelsize::Tuple{Int, Int},
)
    length(drain) == length(indices) == length(surface_flow_width) ||
        throw(DimensionMismatch("drain data must match the land domain"))
    drain = collect(drain)
    false_drains = findall(
        index -> !iszero(drain[index]) && iszero(surface_flow_width[index]),
        eachindex(drain),
    )
    drain[false_drains] .= 0
    land_indices = findall(!iszero, drain)
    drain_indices = indices[land_indices]
    reverse_indices = zeros(Int, modelsize)
    for (drain_idx, cartesian_idx) in enumerate(drain_indices)
        reverse_indices[cartesian_idx] = drain_idx
    end
    return NetworkDrain(; indices = drain_indices, reverse_indices, land_indices)
end
