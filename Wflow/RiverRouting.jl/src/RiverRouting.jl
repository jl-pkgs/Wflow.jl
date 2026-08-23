module RiverRouting

using Dates: TimeType, dayofyear, isleapyear
using DelimitedFiles: readdlm
using EnumX: @enumx
using Graphs: inneighbors, ne, outneighbors, topological_sort_by_dfs
using Parameters: @kwdef, @with_kw
using Polyester: @batch
using StaticArrays: SVector, setindex
using Statistics: mean, quantile!

import RiverGraphs:
    EdgeConnectivity,
    NetworkLand,
    NetworkReservoir,
    NetworkRiver,
    NodesAtEdge,
    EdgesAtNode,
    LDD_PIT,
    PCR_DIR,
    adjacent_edges_at_node,
    adjacent_nodes_at_edge,
    flowgraph

@enumx RoutingType kinematic_wave manning_staggered local_inertial
@enumx GwfConductivityProfileType uniform exponential
@enumx VerticalConductivityProfile exponential exponential_constant layered layered_exponential

abstract type AbstractRoutingConfig end
abstract type AbstractDomain end
abstract type AbstractDomainLand end
abstract type AbstractDomainRiver end
abstract type AbstractLandParameters end
abstract type AbstractRiverParameters end
abstract type AbstractSoilModel end
abstract type AbstractSoilParameters end
abstract type AbstractSnowModel end
abstract type AbstractRoutingClock end

include("core_utils.jl")
include("allocation.jl")
include("interfaces.jl")
include("flow_width.jl")
include("routing.jl")
include("utils.jl")
include("timestepping.jl")
include("subsurface/connectivity.jl")
include("subsurface/groundwater.jl")
include("subsurface/lateral_subsurface_flow.jl")
include("subsurface/initialization.jl")
include("subsurface/subsurface_process.jl")
include("subsurface/boundary_conditions.jl")
include("surface/reservoir.jl")
include("surface/floodplain.jl")
include("surface/surface_flow.jl")
include("surface/surface_kinwave.jl")
include("surface/surface_staggered_scheme.jl")
include("surface/surface_process.jl")

end # module RiverRouting
