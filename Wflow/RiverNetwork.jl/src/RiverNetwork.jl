module RiverNetwork

using Accessors: @reset
using Base.Threads: nthreads
using Graphs:
    DiGraph,
    Graph,
    Graphs,
    SimpleDiGraph,
    add_edge!,
    induced_subgraph,
    inneighbors,
    is_cyclic,
    outneighbors,
    topological_sort_by_dfs

export DIRS,
    LDD_PIT,
    PCR_DIR,
    EdgeConnectivity,
    EdgesAtNode,
    NetworkDrain,
    NetworkLand,
    NetworkReservoir,
    NetworkRiver,
    NodesAtEdge,
    active_indices,
    add_vertex_edge_graph!,
    adjacent_edges_at_node,
    adjacent_nodes_at_edge,
    fillnodata_upstream,
    filter_upstream_nodes,
    flowgraph,
    get_drainage_network,
    graph_from_nodes,
    kinwave_set_subdomains,
    network_subdomains,
    set_pit_ldd,
    stream_order,
    subbasins,
    subbasins_order

const LDD_PIT = 5

"Map from PCRaster LDD value to a CartesianIndex."
const PCR_DIR = [
    CartesianIndex(-1, -1),
    CartesianIndex(0, -1),
    CartesianIndex(1, -1),
    CartesianIndex(-1, 0),
    CartesianIndex(0, 0),
    CartesianIndex(1, 0),
    CartesianIndex(-1, 1),
    CartesianIndex(0, 1),
    CartesianIndex(1, 1),
]

include("river_graph.jl")
using .RiverGraph:
    EdgesAtNode,
    NodesAtEdge,
    add_vertex_edge_graph!,
    adjacent_edges_at_node,
    adjacent_nodes_at_edge

include("subdomains.jl")
include("network.jl")

end # module RiverNetwork
