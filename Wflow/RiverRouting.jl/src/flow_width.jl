"""
    set_effective_flowwidth!(width_x, width_y, domain)

Correct effective floodplain flow widths at cell edges by subtracting river widths. For
D8 diagonal directions, half of the river width is assigned to each adjacent edge.
Reservoir outlet edges are closed. All widths are in metres.
"""
function set_effective_flowwidth!(
    width_x::Vector{Float64},
    width_y::Vector{Float64},
    domain::AbstractDomain,
)::Nothing
    (; local_drain_direction, indices) = domain.river.network
    (; edge_indices, reverse_indices) = domain.land.network
    (; flow_width, reservoir_outlet) = domain.river.parameters
    reverse_indices = reverse_indices[indices]

    graph, local_drain_direction = flowgraph(local_drain_direction, indices, PCR_DIR)
    topological_order = topological_sort_by_dfs(graph)
    n_cells = length(width_x)
    for river_idx in topological_order
        downstream_nodes = outneighbors(graph, river_idx)
        isempty(downstream_nodes) && continue
        downstream_idx = only(downstream_nodes)
        river_width = min(flow_width[river_idx], flow_width[downstream_idx])
        direction = PCR_DIR[local_drain_direction[river_idx]]
        land_idx = reverse_indices[river_idx]
        is_reservoir = reservoir_outlet[river_idx]
        if direction == CartesianIndex(1, 1)
            width_x[land_idx] =
                is_reservoir ? 0.0 : max(width_x[land_idx] - 0.5 * river_width, 0.0)
            width_y[land_idx] =
                is_reservoir ? 0.0 : max(width_y[land_idx] - 0.5 * river_width, 0.0)
        elseif direction == CartesianIndex(-1, -1)
            x_down_idx = edge_indices.ind_x_down[land_idx]
            y_down_idx = edge_indices.ind_y_down[land_idx]
            if x_down_idx <= n_cells
                width_y[x_down_idx] =
                    is_reservoir ? 0.0 : max(width_y[x_down_idx] - 0.5 * river_width, 0.0)
            end
            if y_down_idx <= n_cells
                width_x[y_down_idx] =
                    is_reservoir ? 0.0 : max(width_x[y_down_idx] - 0.5 * river_width, 0.0)
            end
        elseif direction == CartesianIndex(1, 0)
            width_y[land_idx] =
                is_reservoir ? 0.0 : max(width_y[land_idx] - river_width, 0.0)
        elseif direction == CartesianIndex(0, 1)
            width_x[land_idx] =
                is_reservoir ? 0.0 : max(width_x[land_idx] - river_width, 0.0)
        elseif direction == CartesianIndex(-1, 0)
            x_down_idx = edge_indices.ind_x_down[land_idx]
            if x_down_idx <= n_cells
                width_y[x_down_idx] =
                    is_reservoir ? 0.0 : max(width_y[x_down_idx] - river_width, 0.0)
            end
        elseif direction == CartesianIndex(0, -1)
            y_down_idx = edge_indices.ind_y_down[land_idx]
            if y_down_idx <= n_cells
                width_x[y_down_idx] =
                    is_reservoir ? 0.0 : max(width_x[y_down_idx] - river_width, 0.0)
            end
        elseif direction == CartesianIndex(1, -1)
            width_y[land_idx] = max(width_y[land_idx] - 0.5 * river_width, 0.0)
            y_down_idx = edge_indices.ind_y_down[land_idx]
            if y_down_idx <= n_cells
                width_x[y_down_idx] =
                    is_reservoir ? 0.0 : max(width_x[y_down_idx] - 0.5 * river_width, 0.0)
            end
        elseif direction == CartesianIndex(-1, 1)
            x_down_idx = edge_indices.ind_x_down[land_idx]
            if x_down_idx <= n_cells
                width_y[x_down_idx] =
                    is_reservoir ? 0.0 : max(width_y[x_down_idx] - 0.5 * river_width, 0.0)
            end
            width_x[land_idx] =
                is_reservoir ? 0.0 : max(width_x[land_idx] - 0.5 * river_width, 0.0)
        end
    end
    return nothing
end
