@testitem "integration: unconfined transient 1D" begin
    using StaticArrays: SVector
    include("testing_utils.jl")

    nrow = 1
    ncol = 9
    shape = (nrow, ncol)
    conductivity = 0.0023148148148148147
    top = 150.0
    bottom = 0.0
    specific_yield = 0.15
    cellsize = 500.0
    beta = 1.12
    aquifer_length = cellsize * ncol
    gwf_f = 3.0
    conductivity_profile = Wflow.GwfConductivityProfileType.uniform

    # Domain, geometry
    domain = ones(Bool, shape)
    dx = fill(cellsize, ncol)
    dy = fill(cellsize, nrow)
    indices, reverse_indices = Wflow.active_indices(domain, false)
    connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
    ncell = connectivity.ncell
    xc = collect(range(0.0; stop = aquifer_length - cellsize, step = cellsize))

    # constant head on left boundary, 0 at 0
    variables = Wflow.ConstantHeadVariables(; head = [0.0])
    constanthead = Wflow.ConstantHead(; variables, index = [1])

    variables = Wflow.GroundwaterFlowVariables(;
        n = ncell,
        head = initial_head.(xc),
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
        exfiltwater_cumulative = fill(0.0, ncell),
    )
    parameters = Wflow.GroundwaterFlowParameters(;
        hydraulic_conductivity = fill(conductivity, ncell),
        top = fill(top, ncell),
        bottom = fill(bottom, ncell),
        area = fill(cellsize * cellsize, ncell),
        specific_yield = fill(specific_yield, ncell),
        hydraulic_conductivity_scale_parameter = fill(gwf_f, ncell),
    )

    timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.25)
    gwf_model = Wflow.GroundwaterFlowModel(;
        timestepping,
        parameters,
        variables,
        connectivity,
        constanthead,
    )
    domain = Wflow.Domain()

    N = 1
    n = ncell
    water_table_depth = @. 1000.0 * (gwf_model.parameters.top - gwf_model.variables.head)
    soil_model = init_sbm_soil_model(
        n,
        N;
        # Variables
        unsaturated_layer_thickness = SVector.(water_table_depth),
        unsaturated_layer_depth = SVector.(zeros(n)),
        n_unsatlayers = fill(N, n),
        water_table_depth,
        # Parameters
        maximum_number_of_layers = N,
        number_of_layers = fill(1, n),
        theta_s = fill(0.45, n),
        theta_r = fill(0.05, n),
    )

    time = 1.728e6
    t = 0.0
    (; alpha_coefficient) = gwf_model.timestepping
    while t < time
        global t
        gwf_model.variables.q_net .= 0.0
        dt_s = Wflow.stable_timestep(gwf_model, conductivity_profile, alpha_coefficient)
        dt_s = Wflow.check_timestepsize(dt_s, t, time)
        Wflow.update_fluxes!(gwf_model, domain, conductivity_profile, dt_s)
        Wflow.update_head!(gwf_model, soil_model, dt_s)
        t = t + dt_s
        t += dt_s
        # Gradient dh/dx is positive, all flow to the left
        @test all(diff(gwf_model.variables.head) .> 0.0)
    end
end

@testitem "integration: unconfined transient 1D, exponential conductivity" begin
    using StaticArrays: SVector
    include("testing_utils.jl")
    dt = 86400.0
    nrow = 1
    ncol = 9
    shape = (nrow, ncol)
    conductivity = 2.3148148148148148e-6
    top = 150.0
    bottom = 0.0
    specific_yield = 0.15
    cellsize = 500.0
    beta = 1.12
    aquifer_length = cellsize * ncol
    gwf_f = 3.0
    conductivity_profile = Wflow.GwfConductivityProfileType.exponential

    # Domain, geometry
    domain = ones(Bool, shape)
    dx = fill(cellsize, ncol)
    dy = fill(cellsize, nrow)
    indices, reverse_indices = Wflow.active_indices(domain, false)
    connectivity = Wflow.Connectivity(indices, reverse_indices, dx, dy)
    ncell = connectivity.ncell
    xc = collect(range(0.0; stop = aquifer_length - cellsize, step = cellsize))

    # constant head on left boundary, 0 at 0
    variables = Wflow.ConstantHeadVariables(; head = [0.0])
    constanthead = Wflow.ConstantHead(; variables, index = [1])

    variables = Wflow.GroundwaterFlowVariables(;
        n = ncell,
        head = initial_head.(xc),
        conductance = fill(0.0, connectivity.nconnection),
        storage = fill(0.0, ncell),
        q_net = fill(0.0, ncell),
    )
    parameters = Wflow.GroundwaterFlowParameters(;
        hydraulic_conductivity = fill(conductivity, ncell),
        top = fill(top, ncell),
        bottom = fill(bottom, ncell),
        area = fill(cellsize * cellsize, ncell),
        specific_yield = fill(specific_yield, ncell),
        hydraulic_conductivity_scale_parameter = fill(gwf_f, ncell),
    )

    timestepping = Wflow.TimeStepping(; alpha_coefficient = 0.25)
    gwf_model = Wflow.GroundwaterFlowModel(;
        timestepping,
        parameters,
        variables,
        connectivity,
        constanthead,
    )
    domain = Wflow.Domain()

    N = 1
    n = ncell
    water_table_depth = @. 1000.0 * (gwf_model.parameters.top - gwf_model.variables.head)
    soil_model = init_sbm_soil_model(
        n,
        N;
        # Variables
        unsaturated_layer_thickness = SVector.(water_table_depth),
        unsaturated_layer_depth = SVector.(zeros(n)),
        n_unsatlayers = fill(N, n),
        water_table_depth,
        # Parameters
        maximum_number_of_layers = N,
        number_of_layers = fill(1, n),
        theta_s = fill(0.45, n),
        theta_r = fill(0.05, n),
    )

    time = 1.728e6
    t = 0.0
    (; alpha_coefficient) = gwf_model.timestepping
    while t < time
        global t
        gwf_model.variables.q_net .= 0.0
        dt_s = Wflow.stable_timestep(gwf_model, conductivity_profile, alpha_coefficient)
        dt_s = Wflow.check_timestepsize(dt_s, t, time)
        Wflow.update_fluxes!(gwf_model, domain, conductivity_profile, dt_s)
        Wflow.update_head!(gwf_model, soil_model, dt_s)
        t += dt_s
        # Gradient dh/dx is positive, all flow to the left
        @test all(diff(gwf_model.variables.head) .> 0.0)
    end
end
