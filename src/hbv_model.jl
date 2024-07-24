"""
    initialize_hbv_model(config::Config)

Initial part of the SBM model concept. Reads the input settings and data as defined in the
Config object. Will return a Model that is ready to run.
"""
function initialize_hbv_model(config::Config)
    # unpack the paths to the netCDF files
    static_path = input_path(config, config.input.path_static)

    reader = prepare_reader(config)
    clock = Clock(config, reader)
    dt = clock.dt

    do_reservoirs = get(config.model, "reservoirs", false)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    do_pits = get(config.model, "pits", false)::Bool
    set_kquickflow = get(config.model, "set_kquickflow", false)::Bool

    kw_river_tstep = get(config.model, "kw_river_tstep", 0)
    kw_land_tstep = get(config.model, "kw_land_tstep", 0)
    kinwave_it = get(config.model, "kin_wave_iteration", false)::Bool

    nc = NCDataset(static_path)

    subcatch_2d = ncread(nc, config, "subcatchment"; optional=false, allow_missing=true)
    # indices based on catchment
    inds, rev_inds = active_indices(subcatch_2d, missing)
    n = length(inds)

    kw = (; sel=inds, type=Float)
    cfmax = ncread(nc, config, "vertical.cfmax"; kw..., defaults=3.75653) .* (dt / basetimestep)
    tt = ncread(nc, config, "vertical.tt"; kw..., defaults=-1.41934)
    tti = ncread(nc, config, "vertical.tti"; kw..., defaults=1.0)
    ttm = ncread(nc, config, "vertical.ttm"; kw..., defaults=-1.41934)
    whc = ncread(nc, config, "vertical.whc"; kw..., defaults=0.1)

    # glacier parameters
    g_tt = ncread(nc, config, "vertical.g_tt"; kw..., defaults=0.0, fill=0.0)
    g_cfmax = ncread(nc, config, "vertical.g_cfmax"; kw..., defaults=3.0, fill=0.0) .* (dt / basetimestep)
    g_sifrac = ncread(nc, config, "vertical.g_sifrac"; kw..., defaults=0.001, fill=0.0) .* (dt / basetimestep)
    glacierfrac = ncread(nc, config, "vertical.glacierfrac"; kw..., defaults=0.0, fill=0.0)
    glacierstore = ncread(nc, config, "vertical.glacierstore"; kw..., defaults=5500.0, fill=0.0)

    fc = ncread(nc, config, "vertical.fc"; kw..., defaults=260.0)
    betaseepage = ncread(nc, config, "vertical.betaseepage"; kw..., defaults=1.8)
    lp = ncread(nc, config, "vertical.lp"; kw..., defaults=0.53)
    k4 = ncread(nc, config, "vertical.k4"; kw..., defaults=0.02307) .* (dt / basetimestep)
    kquickflow = ncread(nc, config, "vertical.kquickflow"; kw..., defaults=0.09880) .* (dt / basetimestep)
    suz = ncread(nc, config, "vertical.suz"; kw..., defaults=100.0)
    k0 = ncread(nc, config, "vertical.k0"; kw..., defaults=0.30) .* (dt / basetimestep)
    khq = ncread(nc, config, "vertical.khq"; kw..., defaults=0.09880) .* (dt / basetimestep)
    hq = ncread(nc, config, "vertical.hq"; kw..., defaults=3.27) .* (dt / basetimestep)
    alphanl = ncread(nc, config, "vertical.alphanl"; kw..., defaults=1.1)
    perc = ncread(nc, config, "vertical.perc"; kw..., defaults=0.4) .* (dt / basetimestep)
    cfr = ncread(nc, config, "vertical.cfr"; kw..., defaults=0.05)
    pcorr = ncread(nc, config, "vertical.pcorr"; kw..., defaults=1.0)
    rfcf = ncread(nc, config, "vertical.rfcf"; kw..., defaults=1.0)
    sfcf = ncread(nc, config, "vertical.sfcf"; kw..., defaults=1.0)
    cflux = ncread(nc, config, "vertical.cflux"; kw..., defaults=2.0) .* (dt / basetimestep)
    icf = ncread(nc, config, "vertical.icf"; kw..., defaults=2.0)
    cevpf = ncread(nc, config, "vertical.cevpf"; kw..., defaults=1.0)
    epf = ncread(nc, config, "vertical.epf"; kw..., defaults=1.0)
    ecorr = ncread(nc, config, "vertical.ecorr"; kw..., defaults=1.0)

    # read x, y coordinates and calculate cell length [m]
    y_nc = read_y_axis(nc)
    x_nc = read_x_axis(nc)
    y = permutedims(repeat(y_nc, outer=(1, length(x_nc))))[inds]
    cellength = abs(mean(diff(x_nc)))

    sizeinmetres = get(config.model, "sizeinmetres", false)::Bool
    xl, yl = cell_lengths(y, cellength, sizeinmetres)

    threshold = fc .* lp

    hbv = HBV{Float}(;
        dt=Float(tosecond(dt)),
        n, fc, betaseepage,
        lp,
        threshold,
        k4,
        kquickflow=set_kquickflow ? kquickflow :
                   pow.(khq, 1.0 .+ alphanl) .* pow.(hq, -alphanl),
        suz, k0,
        khq, hq,
        alphanl,
        perc, cfr, pcorr,
        rfcf, sfcf, cflux, icf, cevpf, epf, ecorr,
        tti, tt, ttm,
        cfmax,
        whc,

        # glacier parameters
        g_tt,
        g_sifrac,
        g_cfmax,
        glacierstore,
        glacierfrac,

        # default (cold) states:
        interceptionstorage=zeros(Float, n),
        snow=zeros(Float, n),
        snowwater=zeros(Float, n),
        soilmoisture=copy(fc),
        upperzonestorage=0.2 .* fc,
        lowerzonestorage=1.0 ./ (3.0 .* k4),

        # variables:
        precipitation=fill(mv, n),
        temperature=fill(mv, n),
        potential_evaporation=fill(mv, n),
        potsoilevap=fill(mv, n),
        soilevap=fill(mv, n),
        intevap=fill(mv, n),
        actevap=fill(mv, n),
        rainfallplusmelt=fill(mv, n),
        directrunoff=fill(mv, n),
        hbv_seepage=fill(mv, n),
        in_upperzone=fill(mv, n),
        quickflow=fill(mv, n),
        real_quickflow=fill(mv, n),
        percolation=fill(mv, n),
        capflux=fill(mv, n),
        baseflow=fill(mv, n),
        runoff=fill(mv, n),
    )

    modelsize_2d = size(subcatch_2d)
    river_2d = ncread(nc, config, "river_location"; optional=false, type=Bool, fill=false)
    river = river_2d[inds]

    riverwidth_2d = ncread(nc, config, "lateral.river.width"; optional=false, type=Float, fill=0)
    riverwidth = riverwidth_2d[inds]

    riverlength_2d = ncread(nc, config, "lateral.river.length"; optional=false, type=Float, fill=0)
    riverlength = riverlength_2d[inds]

    inds_riv, rev_inds_riv = active_indices(river_2d, 0)
    nriv = length(inds_riv)

    # reservoirs
    pits = zeros(Bool, modelsize_2d)
    if do_reservoirs
        reservoirs, resindex, reservoir, pits =
            initialize_simple_reservoir(config, nc, inds_riv, nriv, pits, tosecond(dt))
    else
        reservoir = ()
        reservoirs = nothing
        resindex = fill(0, nriv)
    end

    # lakes
    if do_lakes
        lakes, lakeindex, lake, pits =
            initialize_lake(config, nc, inds_riv, nriv, pits, tosecond(dt))
    else
        lake = ()
        lakes = nothing
        lakeindex = fill(0, nriv)
    end

    ldd_2d = ncread(nc, config, "ldd"; optional=false, allow_missing=true)
    ldd = ldd_2d[inds]
    if do_pits
        pits_2d = ncread(nc, config, "pits"; optional=false, type=Bool, fill=false)
        ldd = set_pit_ldd(pits_2d, ldd, inds)
    end

    landslope = ncread(nc, config, "lateral.land.slope"; optional=false, kw...)
    clamp!(landslope, 0.00001, Inf)

    dl = map(detdrainlength, ldd, xl, yl)
    dw = (xl .* yl) ./ dl
    olf = initialize_surfaceflow_land(
        nc, config, inds;
        sl=landslope,
        dl=dl,
        width=map(det_surfacewidth, dw, riverwidth, river),
        iterate=kinwave_it,
        tstep=kw_land_tstep,
        dt=dt,
    )

    graph = flowgraph(ldd, inds, pcr_dir)
    riverlength = riverlength_2d[inds_riv]
    riverwidth = riverwidth_2d[inds_riv]
    minimum(riverlength) > 0 || error("river length must be positive on river cells")
    minimum(riverwidth) > 0 || error("river width must be positive on river cells")

    ldd_riv = ldd_2d[inds_riv]
    if do_pits
        ldd_riv = set_pit_ldd(pits_2d, ldd_riv, inds_riv)
    end
    graph_riv = flowgraph(ldd_riv, inds_riv, pcr_dir)

    # the indices of the river cells in the land(+river) cell vector
    index_river = filter(i -> !isequal(river[i], 0), 1:n)
    frac_toriver = fraction_runoff_toriver(graph, ldd, index_river, landslope, n)

    rf = initialize_surfaceflow_river(
        nc, config, inds_riv;
        dl=riverlength,
        width=riverwidth,
        reservoir_index=resindex,
        reservoir=reservoirs,
        lake_index=lakeindex,
        lake=lakes,
        iterate=kinwave_it,
        tstep=kw_river_tstep,
        dt=dt,
    )

    # setup subdomains for the land and river kinematic wave domain, if nthreads = 1
    # subdomain is equal to the complete domain
    toposort = topological_sort_by_dfs(graph)
    toposort_riv = topological_sort_by_dfs(graph_riv)
    index_pit_land = findall(x -> x == 5, ldd)
    index_pit_river = findall(x -> x == 5, ldd_riv)
    streamorder = stream_order(graph, toposort)
    min_streamorder_land = get(config.model, "min_streamorder_land", 5)
    subbas_order, indices_subbas, topo_subbas = kinwave_set_subdomains(
        graph,
        toposort,
        index_pit_land,
        streamorder,
        min_streamorder_land,
    )
    min_streamorder_river = get(config.model, "min_streamorder_river", 6)
    subriv_order, indices_subriv, topo_subriv = kinwave_set_subdomains(
        graph_riv,
        toposort_riv,
        index_pit_river,
        streamorder[index_river],
        min_streamorder_river,
    )

    modelmap = (vertical=hbv, lateral=(land=olf, river=rf))
    indices_reverse = (
        land=rev_inds,
        river=rev_inds_riv,
        reservoir=isempty(reservoir) ? nothing : reservoir.reverse_indices,
        lake=isempty(lake) ? nothing : lake.reverse_indices,
    )
    writer = prepare_writer(config, modelmap, indices_reverse, x_nc, y_nc, nc)
    close(nc)

    # for each domain save:
    # - the directed acyclic graph (graph),
    # - the traversion order (order),
    # - upstream_nodes,
    # - subdomains for the kinematic wave domains for parallel execution (execution order of
    #   subbasins (subdomain_order), traversion order per subbasin (topo_subdomain) and
    #   Vector indices per subbasin matching the traversion order of the complete domain
    #   (indices_subdomain))
    # - the indices that map it back to the two dimensional grid (indices)

    # for the land domain the x and y length [m] of the grid cells are stored
    # for reservoirs and lakes indices information is available from the initialization
    # functions
    land = (
        graph=graph,
        upstream_nodes=filter_upsteam_nodes(graph, pits[inds]),
        subdomain_order=subbas_order,
        topo_subdomain=topo_subbas,
        indices_subdomain=indices_subbas,
        order=toposort,
        indices=inds,
        reverse_indices=rev_inds,
        area=xl .* yl,
    )
    river = (
        graph=graph_riv,
        upstream_nodes=filter_upsteam_nodes(graph_riv, pits[inds_riv]),
        subdomain_order=subriv_order,
        topo_subdomain=topo_subriv,
        indices_subdomain=indices_subriv,
        order=toposort_riv,
        indices=inds_riv,
        reverse_indices=rev_inds_riv,
    )

    model = Model(
        config,
        (; land, river, reservoir, lake, index_river, frac_toriver),
        (subsurface=nothing, land=olf, river=rf),
        hbv,
        clock,
        reader,
        writer,
        HbvModel(),
    )

    model = set_states(model)

    return model
end

function update(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:HbvModel}
    @unpack lateral, vertical, network, clock, config = model

    inds_riv = network.index_river

    # vertical hbv concept is updated until snow state, after that (optional)
    # snow transport is possible
    update_until_snow(vertical, config)

    # lateral snow transport
    if get(config.model, "masswasting", false)::Bool
        lateral_snow_transport!(
            vertical.snow,
            vertical.snowwater,
            lateral.land.sl,
            network.land,
        )
    end

    # update vertical hbv concept
    update_after_snow(vertical, config)

    surface_routing(model)

    return model
end

function set_states(model::Model{N,L,V,R,W,T}) where {N,L,V,R,W,T<:HbvModel}
    @unpack lateral, config = model

    reinit = get(config.model, "reinit", true)::Bool
    do_lakes = get(config.model, "lakes", false)::Bool
    # read and set states in model object if reinit=true
    if reinit == false
        instate_path = input_path(config, config.state.path_input)
        @info "Set initial conditions from state file `$instate_path`."
        set_states(instate_path, model; type=Float)
        # update kinematic wave volume for river and land domain
        @unpack lateral = model
        # makes sure land cells with zero flow width are set to zero q and h
        for i in eachindex(lateral.land.width)
            if lateral.land.width[i] <= 0.0
                lateral.land.q[i] = 0.0
                lateral.land.h[i] = 0.0
            end
        end
        lateral.land.volume .= lateral.land.h .* lateral.land.width .* lateral.land.dl
        lateral.river.volume .= lateral.river.h .* lateral.river.width .* lateral.river.dl

        if do_lakes
            # storage must be re-initialized after loading the state with the current
            # waterlevel otherwise the storage will be based on the initial water level
            lakes = lateral.river.lake
            lakes.storage .=
                initialize_storage(lakes.storfunc, lakes.area, lakes.waterlevel, lakes.sh)
        end
    else
        @info "Set initial conditions from default values."
    end

    return model
end
