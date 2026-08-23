@testitem "kinematic wave overland flow" begin
    using NCDatasets: NCDataset
    using Graphs: topological_sort_by_dfs

    dt_sec = 86400.0
    ldd_MISSING_VALUE = 255

    # read the staticmaps into memory
    nc = NCDataset(normpath(@__DIR__, "data/input/staticmaps-rhine.nc"))
    # helper function to get the axis order and directionality right
    read_right(nc, var) = reverse(permutedims(Array(nc[var])); dims = 2)
    ldd_2d = read_right(nc, "ldd")

    inds, _ = Wflow.active_indices(ldd_2d, ldd_MISSING_VALUE)
    n = length(inds)

    # take out only the active cells
    ldd = ldd_2d[inds]
    slope = read_right(nc, "slope")[inds]
    N = read_right(nc, "N")[inds]
    Qold = read_right(nc, "Qold")[inds]
    Bw = read_right(nc, "Bw")[inds]
    waterlevel = read_right(nc, "waterlevel")[inds]
    DCL = read_right(nc, "DCL")[inds]
    close(nc)

    # create the directed acyclic graph from the drainage direction array
    graph, ldd = Wflow.flowgraph(ldd, inds, Wflow.PCR_DIR)
    # a topological sort is used for visiting nodes in order from upstream to downstream
    toposort = topological_sort_by_dfs(graph)
    sink = toposort[end]
    @test ldd[sink] == Wflow.LDD_PIT  # the most downstream node must be a sink

    # calculate parameters of kinematic wave
    q = 0.000001
    beta = Wflow.BETA_KINWAVE
    AlpPow = (2.0 / 3.0) * beta
    AlpTermR = (N ./ sqrt.(slope)) .^ beta
    P = Bw + (2.0 * waterlevel)
    alpha = AlpTermR .* P .^ AlpPow

    Q = zeros(n)
    Q = Wflow.kin_wave!(Q, graph, toposort, Qold, q, alpha, DCL, dt_sec)

    @test sum(Q) ≈ 2.957806043289641e6
    @test Q[toposort[1]] ≈ 0.007260052312634069
    @test Q[toposort[n - 100]] ≈ 3945.762718338739
    @test Q[sink] ≈ 4131.101474418251
end

@testitem "unit: kinematic_wave_ssf" begin
    using StaticArrays: SVector
    include("testing_utils.jl")

    ### Shared values
    n = 1
    N = 4
    actual_layer_thickness = [SVector(0.1, 0.3, 0.8, 0.8)]
    cumulative_layer_depth = [SVector(0.0, 0.1, 0.4, 1.2, 2.0)]
    maximum_number_of_layers = 4
    number_of_layers = [4]
    theta_s = [0.48642662167549133]
    theta_r = [0.11939866840839386]
    theta_fc = [0.28219206182657536]
    ssfin = 0.0

    ssf_prev = 0.30038365579798126
    zi_prev = 0.0005198340870375973
    q_net = 0.005618827458801466
    slope = 0.4522336721420288
    sy = 0.20423455984891598
    d = 2.0
    dt = 86400.0
    dx = 1117.0150713112287
    dw = 517.495693771673
    ssfmax = 0.0009215296489248933
    kh_profile = Wflow.KhExponential([0.002379589787235966], [1.0141291422769427])
    i = 1

    soil_model = init_sbm_soil_model(
        n,
        N;
        # Variables
        unsaturated_layer_thickness = [SVector(0.1, 0.3, 0.11983408703759733, NaN)],
        unsaturated_layer_depth = [
            SVector(0.0001909439890049523, 0.01627933934181815, 0.019508197676020186, 0.0),
        ],
        n_unsatlayers = [3],
        water_table_depth = [0.5198340870375974],
        # Parameters
        maximum_number_of_layers,
        cumulative_layer_depth,
        number_of_layers,
        theta_s,
        theta_r,
        theta_fc,
        actual_layer_thickness,
    )

    # Case: !(ssfin + ssf_prev ≈ 0.0 && qnet <= 0)
    # Case: !(zi > d)
    ssf, water_table_depth, exfilt, net_flux = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        q_net,
        slope,
        sy,
        d,
        dt,
        dx,
        dw,
        ssfmax,
        kh_profile,
        soil_model,
        i,
    )
    @test ssf ≈ 0.23130576097772237
    @test water_table_depth ≈ 0.1656875455413981
    @test exfilt ≈ 0.0
    @test net_flux ≈ -3.904277181728481e-7

    # Case: ssfin + ssf_prev ≈ 0.0 && q_net <= 0
    ssf_prev = 0.0
    q_net = 0.0
    zi_prev = 0.0
    ssfmax = 0.0009215296489248933
    kh_profile = Wflow.KhExponential([0.002379589787235966], [1.0141291422769427])
    ssf, water_table_depth, exfilt, sy_d = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        q_net,
        slope,
        sy,
        d,
        dt,
        dx,
        dw,
        ssfmax,
        kh_profile,
        soil_model,
        i,
    )
    @test iszero(ssf)
    @test water_table_depth == d
    @test iszero(exfilt)
    @test sy_d ≈ 0.0

    soil = init_sbm_soil_model(
        n,
        N;
        # Variables
        unsaturated_layer_thickness = [SVector(0.1, 0.3, 0.348312461531486, NaN)],
        unsaturated_layer_depth = [
            SVector(0.0001909439890049523, 0.01627933934181815, 0.058425012193036086, 0.0),
        ],
        n_unsatlayers = [3],
        water_table_depth = [0.748312461531486],
        # Parameters
        maximum_number_of_layers,
        cumulative_layer_depth,
        number_of_layers,
        theta_s,
        theta_r,
        theta_fc,
        actual_layer_thickness,
    )

    ssf_prev = 0.627032986563781
    zi_prev = 0.748312461531486
    q_net = 0.008957349820205272
    slope = 0.4522336721420288
    sy = 0.20423455984891598
    d = 2.0
    dt = 86400.0
    dx = 1117.0150713112287
    dw = 517.495693771673
    ssfmax = 0.0017762382461437296
    kh_profile = Wflow.KhExponentialConstant(kh_profile, [0.2])
    i = 1

    ssf, water_table_depth, exfilt, net_flux = Wflow.kinematic_wave_ssf(
        ssfin,
        ssf_prev,
        zi_prev,
        q_net,
        slope,
        sy,
        d,
        dt,
        dx,
        dw,
        ssfmax,
        kh_profile,
        soil,
        i,
    )

    @test ssf ≈ 0.5171363105669935
    @test water_table_depth ≈ 1.1202203724020348
    @test exfilt ≈ 0.0
    @test net_flux ≈ -8.791255611224121e-7
end
