"Read a routing input through an application-provided data backend."
function ncread end

"Resolve an application input path."
function input_path end

"Return whether water-demand allocation is enabled."
function do_water_demand end

"Update the soil water table and return water-table change and exfiltration."
function water_table_change end

"Update unsaturated soil layers after a water-table change."
function update_ustorelayerdepth! end

"Initialize a river allocation component for an application configuration."
river_allocation(config, n::Int) =
    do_water_demand(config) ? AllocationRiverModel(; n) : NoAllocationRiverModel(n)
