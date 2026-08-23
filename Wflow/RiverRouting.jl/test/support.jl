using Test
using Parameters: @kwdef, @with_kw
using RiverRouting
using RiverGraphs

@with_kw struct TestLandParameters <: RiverRouting.AbstractLandParameters
    x_length::Vector{Float64} = Float64[]
    y_length::Vector{Float64} = Float64[]
    area::Vector{Float64} = Float64[]
    flow_width::Vector{Float64} = Float64[]
    surface_flow_width::Vector{Float64} = Float64[]
    flow_length::Vector{Float64} = Float64[]
    flow_fraction_to_river::Vector{Float64} = Float64[]
    slope::Vector{Float64} = Float64[]
    reservoir_outlet::Vector{Bool} = Bool[]
    reservoir_coverage::Vector{Bool} = Bool[]
    river_location::Vector{Bool} = Bool[]
    river_fraction::Vector{Float64} = Float64[]
    water_fraction::Vector{Float64} = Float64[]
end

@with_kw struct TestRiverParameters <: RiverRouting.AbstractRiverParameters
    flow_width::Vector{Float64} = Float64[]
    flow_length::Vector{Float64} = Float64[]
    slope::Vector{Float64} = Float64[]
    reservoir_outlet::Vector{Bool} = Bool[]
    reservoir_coverage::Vector{Bool} = Bool[]
    cell_area::Vector{Float64} = Float64[]
end

@kwdef struct TestDomainLand <: RiverRouting.AbstractDomainLand
    network::RiverGraphs.NetworkLand = RiverGraphs.NetworkLand()
    parameters::TestLandParameters = TestLandParameters()
end

@kwdef struct TestDomainRiver <: RiverRouting.AbstractDomainRiver
    network::RiverGraphs.NetworkRiver = RiverGraphs.NetworkRiver()
    parameters::TestRiverParameters = TestRiverParameters()
end

@kwdef struct TestDomainReservoir
    network::RiverGraphs.NetworkReservoir = RiverGraphs.NetworkReservoir()
end

@kwdef struct TestDomainDrain
    network::RiverGraphs.NetworkDrain = RiverGraphs.NetworkDrain()
end

@kwdef struct TestDomain <: RiverRouting.AbstractDomain
    land::TestDomainLand = TestDomainLand()
    river::TestDomainRiver = TestDomainRiver()
    reservoir::TestDomainReservoir = TestDomainReservoir()
    drain::TestDomainDrain = TestDomainDrain()
end

function homogeneous_aquifer(nrow::Int, ncol::Int)
    active = ones(Bool, nrow, ncol)
    indices, reverse_indices = RiverGraphs.active_indices(active, false)
    connectivity = RiverRouting.Connectivity(
        indices,
        reverse_indices,
        fill(10.0, ncol),
        fill(10.0, nrow),
    )
    n_cells = connectivity.ncell
    constant_head = RiverRouting.ConstantHead(;
        variables = RiverRouting.ConstantHeadVariables(; head = Float64[]),
        index = Int[],
    )
    parameters = RiverRouting.GroundwaterFlowParameters(;
        hydraulic_conductivity = fill(10.0 / 86400.0, n_cells),
        top = fill(10.0, n_cells),
        bottom = zeros(n_cells),
        area = fill(100.0, n_cells),
        specific_yield = fill(0.15, n_cells),
        hydraulic_conductivity_scale_parameter = fill(3.0, n_cells),
    )
    variables = RiverRouting.GroundwaterFlowVariables(;
        n = n_cells,
        head = [0.0, 7.5, 20.0],
        conductance = zeros(connectivity.nconnection),
        storage = zeros(n_cells),
        q_net = zeros(n_cells),
        exfiltwater_cumulative = zeros(n_cells),
    )
    return RiverRouting.GroundwaterFlowModel(;
        timestepping = RiverRouting.TimeStepping(),
        parameters,
        variables,
        connectivity,
        constanthead = constant_head,
    )
end
