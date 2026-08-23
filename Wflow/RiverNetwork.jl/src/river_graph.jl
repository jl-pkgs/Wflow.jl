module RiverGraph

using Graphs: SimpleDiGraph, add_edge!, add_vertex!, dst, edges, nv, src, vertices

"Struct for storing source `src` node and destination `dst` node of an edge."
@kwdef struct NodesAtEdge
    src::Vector{Int} = Int[]
    dst::Vector{Int} = Int[]
end

"Struct for storing source `src` edge and destination `dst` edge of a node."
@kwdef struct EdgesAtNode
    src::Vector{Vector{Int}} = Vector{Int}[]
    dst::Vector{Vector{Int}} = Vector{Int}[]
end

"Initialize `NodesAtEdge` and add outlet edges for `pit_nodes` to `graph`."
function NodesAtEdge(graph::SimpleDiGraph{Int}, pit_nodes::Vector{Int})
    add_vertex_edge_graph!(graph, pit_nodes)
    return NodesAtEdge(; adjacent_nodes_at_edge(graph)...)
end

"Initialize `EdgesAtNode` from a directed `graph` and its `nodes_at_edge`."
function EdgesAtNode(graph::SimpleDiGraph{Int}, nodes_at_edge::NodesAtEdge)
    return EdgesAtNode(; adjacent_edges_at_node(graph, nodes_at_edge)...)
end

"""
    adjacent_nodes_at_edge(graph)

Return the source node `src` and destination node `dst` of each edge of a directed `graph`.
"""
function adjacent_nodes_at_edge(
    graph::SimpleDiGraph{Int},
)::NamedTuple{(:src, :dst), Tuple{Vector{Int}, Vector{Int}}}
    graph_edges = collect(edges(graph))
    return (src = src.(graph_edges), dst = dst.(graph_edges))
end

"""
    adjacent_edges_at_node(graph, nodes_at_edge)

Return the source edge `src` and destination edge `dst` of each node of a directed `graph`.
"""
function adjacent_edges_at_node(
    graph::SimpleDiGraph{Int},
    nodes_at_edge,
)::NamedTuple{(:src, :dst), Tuple{Vector{Vector{Int}}, Vector{Vector{Int}}}}
    nodes = vertices(graph)
    source_edges = Vector{Int}[]
    destination_edges = Vector{Int}[]
    for node_idx in 1:nv(graph)
        push!(source_edges, findall(isequal(nodes[node_idx]), nodes_at_edge.dst))
        push!(destination_edges, findall(isequal(nodes[node_idx]), nodes_at_edge.src))
    end
    return (src = source_edges, dst = destination_edges)
end

"Add `vertex` and `edge` to `pits` of a directed `graph`."
function add_vertex_edge_graph!(graph::SimpleDiGraph{Int}, pits::Vector{Int})::Nothing
    n_nodes = nv(graph)
    for (pit_idx, node_idx) in enumerate(pits)
        add_vertex!(graph)
        add_edge!(graph, node_idx, n_nodes + pit_idx)
    end
    return nothing
end

end # module RiverGraph
