using Graphs: DiGraph, Graph, add_edge!, has_edge, nv, topological_sort_by_dfs
using RiverNetwork
using Test

@testset "RiverNetwork" begin
    @testset "network construction" begin
        subcatchment = fill(1, 1, 4)
        drainage = UInt8[8 8 8 5]
        land = NetworkLand(subcatchment, drainage; nodata = missing)

        @test land.indices == CartesianIndex.(Ref(1), 1:4)
        @test land.reverse_indices == reshape(1:4, 1, 4)
        @test land.order == [1, 2, 3, 4]
        @test land.streamorder == ones(Int, 4)
        @test has_edge(land.graph, 1, 2)
        @test has_edge(land.graph, 3, 4)

        river = NetworkRiver(Bool[0 1 1 1], drainage, land)
        @test river.land_indices == [2, 3, 4]
        @test river.indices == CartesianIndex.(Ref(1), 2:4)
        @test river.order == [1, 2, 3]

        edge_indices = EdgeConnectivity(land)
        @test edge_indices.ind_y_down == [5, 1, 2, 3]
        @test edge_indices.ind_y_up == [2, 3, 4, 5]

        graph = copy(river.graph)
        nodes_at_edge = NodesAtEdge(graph, [3])
        edges_at_node = EdgesAtNode(graph, nodes_at_edge)
        @test nodes_at_edge.src == [1, 2, 3]
        @test nodes_at_edge.dst == [2, 3, 4]
        @test edges_at_node.src == [Int[], [1], [2], [3]]
        @test edges_at_node.dst == [[1], [2], [3], Int[]]

        reservoir, reservoir_map = NetworkReservoir([0, 0, 1], [0 0 1 1], river)
        @test reservoir_map == [0, 0, 1]
        @test reservoir.indices_outlet == [CartesianIndex(1, 4)]
        @test reservoir.indices_coverage == [[CartesianIndex(1, 3), CartesianIndex(1, 4)]]

        drain = NetworkDrain([1, 1, 0, 1], land.indices, [0.0, 1.0, 1.0, 1.0], (1, 4))
        @test drain.land_indices == [2, 4]
        @test drain.indices == [CartesianIndex(1, 2), CartesianIndex(1, 4)]
    end

    @testset "drainage graph validation" begin
        indices = CartesianIndex.(Ref(1), 1:2)
        @test_throws ErrorException flowgraph(UInt8[8, 2], indices)

        ldd = UInt8[2, 5]
        @test_logs (:warn,) flowgraph(ldd, indices)
        @test ldd == UInt8[5, 5]
    end

    @testset "stream order and subbasins" begin
        graph = DiGraph(16)
        for (source, destination) in (
            (1, 3),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (7, 9),
            (8, 9),
            (9, 4),
            (10, 12),
            (11, 12),
            (12, 16),
            (13, 15),
            (14, 15),
            (15, 16),
            (16, 5),
        )
            add_edge!(graph, source, destination)
        end

        topological_order = topological_sort_by_dfs(graph)
        order = stream_order(graph, topological_order)
        subbasin = subbasins(graph, order, topological_order, 2)
        filled = fillnodata_upstream(graph, topological_order, subbasin, 0)
        subbasin_graph = graph_from_nodes(graph, subbasin, filled)

        @test order == [1, 1, 2, 3, 4, 4, 1, 1, 2, 1, 1, 2, 1, 1, 2, 3]
        @test subbasin == [0, 0, 5, 6, 0, 7, 0, 0, 4, 0, 0, 2, 0, 0, 1, 3]
        @test filled == [5, 5, 5, 6, 7, 7, 4, 4, 4, 2, 2, 2, 1, 1, 1, 3]
        @test nv(subbasin_graph) == 7
    end
end
