"Common supertype for land and river water-allocation components."
abstract type AbstractAllocationModel end

"Disabled river-allocation component."
struct NoAllocationRiverModel <: AbstractAllocationModel
    n::Int
end

"Variables used to couple river routing and water allocation."
@with_kw struct AllocationRiverVariables
    n::Int
    actual_surfacewater_abstraction::Vector{Float64} = zeros(n)
    actual_surfacewater_abstraction_volume::Vector{Float64} = zeros(n)
    available_surfacewater::Vector{Float64} = zeros(n)
    non_irrigation_returnflow::Vector{Float64} = zeros(n)
end

"River water-allocation component."
@with_kw struct AllocationRiverModel <: AbstractAllocationModel
    n::Int
    variables::AllocationRiverVariables = AllocationRiverVariables(; n)
end

get_nonirrigation_returnflow(allocation_model::AllocationRiverModel) =
    allocation_model.variables.non_irrigation_returnflow
get_nonirrigation_returnflow(allocation_model::NoAllocationRiverModel) =
    zeros(allocation_model.n)
