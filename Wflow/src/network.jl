"Initialize the land network from Wflow input data."
function NetworkLand(dataset::NCDataset, config::Config)
    subcatchment = ncread(dataset, config, "subbasin_location__count", Domain)
    active_ids = config.input.subbasin_active_location__count
    if !isnothing(active_ids)
        @info "Only subcatchments with IDs `$(sort!(collect(Set(active_ids))))` are active."
    end
    ldd = ncread(dataset, config, "basin__local_drain_direction", Domain)
    pits =
        config.model.pit__flag ?
        ncread(dataset, config, "basin_pit_location__mask", Domain) : nothing
    return RiverNetwork.NetworkLand(subcatchment, ldd; active_ids, pits_2d = pits)
end

"Set land-network subdomains from Wflow configuration."
network_subdomains(config::Config, network::NetworkLand) =
    network_subdomains(network, config.model.land_streamorder__min_count)

"Read a drainage network from Wflow input data."
function get_drainage_network(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}};
    do_pits::Bool = false,
    logging::Bool = true,
)
    ldd = ncread(dataset, config, "basin__local_drain_direction", Domain; logging)
    pits = do_pits ? ncread(dataset, config, "basin_pit_location__mask", Domain) : nothing
    return RiverNetwork.get_drainage_network(ldd, indices; pits_2d = pits)
end

"Initialize the river network from Wflow input data."
function NetworkRiver(
    dataset::NCDataset,
    config::Config,
    network::NetworkLand;
    do_pits::Bool = false,
)
    river_location = ncread(
        dataset,
        config,
        "river_location__mask",
        Domain;
        sel = network.indices,
        logging = false,
    )
    ldd = ncread(dataset, config, "basin__local_drain_direction", Domain; logging = false)
    pits = do_pits ? ncread(dataset, config, "basin_pit_location__mask", Domain) : nothing
    return RiverNetwork.NetworkRiver(river_location, ldd, network; pits_2d = pits)
end

"Set river-network subdomains from Wflow configuration."
network_subdomains(config::Config, network::NetworkRiver) =
    network_subdomains(network, config.model.river_streamorder__min_count)

"Initialize the reservoir network from Wflow input data."
function NetworkReservoir(dataset::NCDataset, config::Config, network::NetworkRiver)
    outlet_ids = ncread(
        dataset,
        config,
        "reservoir_location__count",
        Routing;
        sel = network.indices,
        logging = false,
    )
    coverage = ncread(dataset, config, "reservoir_area__count", Routing; logging = false)
    return RiverNetwork.NetworkReservoir(outlet_ids, coverage, network)
end

"Initialize the groundwater-drain network from Wflow input data."
function NetworkDrain(
    dataset::NCDataset,
    config::Config,
    indices::Vector{CartesianIndex{2}},
    surface_flow_width::Vector{Float64},
    modelsize::Tuple{Int, Int},
)
    drain = ncread(dataset, config, "land_drain_location__mask", Routing; sel = indices)
    network = RiverNetwork.NetworkDrain(drain, indices, surface_flow_width, modelsize)
    n_removed = count(!iszero, drain) - length(network.indices)
    n_removed > 0 &&
        @info "$n_removed drain locations are removed that occur where overland flow is not possible (overland flow width is zero)"
    return network
end
