"""
    scurve(x, a, b, c)

Sigmoid "S"-shaped curve.

# Arguments
- `x::Real`: input
- `a::Real`: determines the center level
- `b::Real`: determines the amplitude of the curve (range: (0, b⁻¹))
- `c::Real`: determines the steepness or "stepwiseness" of the curve.
             The higher c the sharper the function. A negative c reverses the function.
"""
function scurve(x::Real, a::Real, b::Real, c::Real)::Real
    s = inv(b + exp(-c * (x - a)))
    return s
end

function to_enumx(T, i::Int)
    options = instances(T)
    n_options = length(options)
    if 1 ≤ i ≤ n_options
        return options[i]
    else
        options_repr = repr(MIME("text/plain"), T)
        throw(
            error(
                "Cannot convert $i to $T, there are only $n_options options:\n$options_repr.",
            ),
        )
    end
end

"Get active indices of `key` (standard name or model path) by prefix string matching."
function active_indices(domain::Domain, key::AbstractString)::Vector{CartesianIndex{2}}
    if startswith(key, domain_parameter_map["reservoir"])
        return domain.reservoir.network.indices_outlet
    elseif startswith(key, domain_parameter_map["river"])
        return domain.river.network.indices
    elseif startswith(key, domain_parameter_map["drain"])
        return domain.drain.network.indices
    else
        return domain.land.network.indices
    end
end

function lattometres(lat::Real)::Tuple{Float64, Float64}
    m1 = 111132.92     # latitude calculation term 1
    m2 = -559.82       # latitude calculation term 2
    m3 = 1.175         # latitude calculation term 3
    m4 = -0.0023       # latitude calculation term 4
    p1 = 111412.84     # longitude calculation term 1
    p2 = -93.5         # longitude calculation term 2
    p3 = 0.118         # longitude calculation term 3

    # Calculate the length of a degree of latitude and longitude in meters
    latlen = m1 + (m2 * cosd(2.0 * lat)) + (m3 * cosd(4.0 * lat)) + (m4 * cosd(6.0 * lat))
    longlen = (p1 * cosd(lat)) + (p2 * cosd(3.0 * lat)) + (p3 * cosd(5.0 * lat))

    return longlen, latlen
end

function cell_lengths(
    y::AbstractVector{<:Real},
    celllength::Real,
    cell_length_in_meter::Bool,
)::Tuple{Vector{Float64}, Vector{Float64}}
    n = length(y)
    xl = fill(MISSING_VALUE, n)
    yl = fill(MISSING_VALUE, n)
    if cell_length_in_meter
        xl .= celllength
        yl .= celllength
    else
        for i in 1:n
            longlen, latlen = lattometres(y[i])
            xl[i] = longlen * celllength
            yl[i] = latlen * celllength
        end
    end
    return xl, yl
end

"""
    set_states!(instate_path, model, state_ncnames; <keyword arguments>)

Read states contained in `Dict` `state_ncnames` from netCDF file located in `instate_path`,
and set states in `model` object. Active cells are selected with the corresponding network's
(`Vector{CartesianIndex}`) from the netCDF file.

# Arguments
- `type = nothing`: type to convert data to after reading. By default no conversion is done.
"""
function set_states!(instate_path::AbstractString, model; dimname = nothing)::Nothing
    (; domain, config, clock, land) = model
    dt_val = tosecond(clock.dt)

    # Check if required states are covered
    state_ncnames = check_states(config)

    # states in netCDF include dim time (one value) at index 3 or 4, 3 or 4 dims are allowed
    NCDataset(instate_path) do ds
        for (state, ncname) in state_ncnames
            metadata = get_metadata(state; model)
            (; unit) = metadata
            @info "Setting initial state from netCDF." ncpath = instate_path ncvarname =
                ncname state unit
            sel = active_indices(domain, state)
            n = length(sel)
            dims = length(dimnames(ds[ncname]))
            # 4 dims, for example (x,y,layer,time) where dim layer is an SVector for soil layers
            if dims == 4
                if dimname == :layer
                    dimensions = (x = :, y = :, layer = :, time = 1)
                else
                    error("Unrecognized dimension name $dimname")
                end
                A = read_standardized(ds, ncname, dimensions)
                A = permutedims(A[sel, :])
                # note that this array is allowed to have missing, since not every land
                # column is `maximum_number_of_layers` layers deep
                if dimname == :layer
                    A = replace!(A, missing => NaN)
                end
                A = apply_unit_and_type_transform!(A, metadata; dt_val)
                # set state in model object
                metadata.lens(model) .= svectorscopy(A, Val{size(A)[1]}())
                # 3 dims (x,y,time)
            elseif dims == 3
                A = read_standardized(ds, ncname, (x = :, y = :, time = 1))
                A = A[sel]
                A = nomissing(A)
                A = apply_unit_and_type_transform!(A, metadata; dt_val)
                # set state in model object, only set active cells ([1:n]) (ignore boundary conditions/ghost points)
                lens = get_metadata(state, typeof(land), Routing; model).lens
                lens(model)[1:n] .= A
            else
                error(
                    "Number of state dims should be 3 or 4, number of dims = ",
                    string(dims),
                )
            end
        end
    end
    return nothing
end

function get_var(config::Config, parameter::AbstractString; optional = true)
    if hasfield(InputSection, Symbol(parameter))
        var = getfield(config.input, Symbol(parameter))
    elseif haskey(config.input._location_maps, parameter)
        var = config.input._location_maps[parameter]
    elseif haskey(config.input.static, parameter)
        var = config.input.static[parameter]
    elseif haskey(config.input.cyclic, parameter)
        var = config.input.cyclic[parameter]
    elseif optional
        var = nothing
    else
        error(
            "Required input model parameter with standard name '$parameter' not set in TOML file",
        )
    end
    return var
end

"""
Apply the affine transform in `var` to the incoming array `A` in place element-wise.
The affine transform consists of a scaling by `scale` and a translation by `offset`.
These operations are only applied when non-trivial.
"""
function apply_affine_transform!(A::AbstractArray, var::InputEntry)
    (; _do_scaling, _scale_scalar, scale, _do_offsetting, _offset_scalar, offset) = var
    if _do_scaling
        if _scale_scalar
            A .*= only(scale)
        else
            A .*= scale
        end
    end
    if _do_offsetting
        if _offset_scalar
            A .+= only(offset)
        else
            A .+= offset
        end
    end
    return A
end

"""
    ncread(nc, config::Config, parameter::AbstractString, model_type; <keyword arguments>)

Read a netCDF variable `var` from file `nc`, based on `config` (parsed TOML file) and the
model `parameter` (standard name) specified in the TOML configuration file. Supports various
keyword arguments to get selections of data in desired types, with or without missing
values.

# Arguments
- `model_type`: The model type (e.g., LandHydrologySBM, SoilLoss, Domain, Routing) used to
        determine the appropriate standard name mapping.
- `sel=nothing`: A selection of indices, such as a `Vector{CartesianIndex}` of active cells,
        to return from the netCDF. By default all cells are returned.
- `logging=true`: Generate a logging message when reading a netCDF variable.
- `metadata`: The metadata of the read parameter or state, obtained from the parameter name by default.
"""
function ncread(
    nc,
    config::Config,
    parameter::AbstractString,
    model_type;
    sel = nothing,
    logging = true,
    metadata = get_metadata(parameter, model_type),
)
    (; default, fill, type, allow_missing, dimname) = metadata
    var = get_var(config, parameter; optional = !isnothing(default))
    dt_val = config.time.timestepsecs
    (; unit) = metadata

    # for optional parameters default values are used.
    if isnothing(var)
        @info "Set `$parameter [$unit]` using default value `$default $unit`."
        @assert !isnothing(default) "Default value required but not available for $parameter (if you see this as a user please open an issue)."
        default = unit_and_type_transform(default, metadata; dt_val)
        if isnothing(dimname)
            return Base.fill(default, length(sel))
        else
            return Base.fill(default, (nc.dim[String(dimname)], length(sel)))
        end
    end

    # dim `time` is also included in `dim_sel`: this allows for cyclic parameters (read
    # first timestep), that is later updated with the `update_cyclic!` function.
    if isnothing(dimname)
        dim_sel = (x = :, y = :, time = 1)
    elseif dimname == :layer
        dim_sel = (x = :, y = :, layer = :, time = 1)
    elseif dimname == :flood_depth
        dim_sel = (x = :, y = :, flood_depth = :, time = 1)
    else
        error("Unrecognized dimension name $dimname")
    end

    if var isa Number
        var = InputEntry(; value = var)
    elseif var isa String
        var = InputEntry(; external_name = var)
    else
        @assert var isa InputEntry
    end

    (; value, layer, scale, offset) = var
    variable_info(var)

    if !isnothing(value)
        @info "Set `$parameter [$unit]` using uniform value `$value $unit` from TOML file."
        A = if isnothing(dimname)
            # set to one uniform value
            Base.fill(only(value), length(sel))
        elseif length(value) == 1
            # set to one uniform value (parameter with third dimension of size 1)
            Base.fill(only(value), (nc.dim[String(dimname)], length(sel)))
        elseif length(value) > 1
            # set to multiple uniform values (parameter with third dimension of size > 1)
            @assert length(value) == nc.dim[String(dimname)]
            repeat(value, 1, length(sel))
        end
        return apply_unit_and_type_transform!(A, metadata; dt_val)
    else
        if logging
            @info "Set `$parameter [$unit]` using netCDF variable `$var`."
        end
        A = read_standardized(nc, variable_name(var), dim_sel)
        if !isnothing(layer)
            # the modifier index is only set in combination with scale and offset for SVectors,
            # provided through the TOML file.
            # if index, scale and offset is provided in the TOML as a list.
            for i in eachindex(layer)
                A[:, :, layer[i]] = A[:, :, layer[i]] .* scale[i] .+ offset[i]
            end
        else
            apply_affine_transform!(A, var)
        end
    end

    # Take out only the active cells
    if !isnothing(sel)
        if isnothing(dimname)
            A = A[sel]
        else
            A = permutedims(A[sel, :])
        end
    end

    if allow_missing
        # Convert to desired type if needed
        A = map(x -> ismissing(x) ? x : type(x), A)
    else
        if isnothing(fill)
            # errors if missing are found
            A = nomissing(A)
            if any(isnan, A)
                error("NaN not allowed in $var")
            end
        else
            # replaces missing with a fill value
            A = nomissing(A, fill)
            # replace also NaN values with the fill value
            replace!(x -> isnan(x) ? fill : x, A)
        end
    end
    return apply_unit_and_type_transform!(A, metadata; dt_val)
end

"""
    get_flow_length(ldd, x_length, y_length)

Return the flow length for a non square grid. Input `ldd` (drainage network), `x_length`
(length of cells in x direction), `y_length` (length of cells in y direction). Output is
flow length.
"""
function get_flow_length(ldd::UInt8, x_length::Real, y_length::Real)::Real
    # take into account non-square cells
    # if ldd is 8 or 2 use y_length
    # if ldd is 4 or 6 use x_length
    if ldd == 2 || ldd == 8
        y_length
    elseif ldd == 4 || ldd == 6
        x_length
    else
        hypot(x_length, y_length)
    end
end

"""
    get_flow_width(ldd, x_length, y_length)

Return the flow width for a non square grid. Input `ldd` (drainage network), `x_length`
(length of cells in x direction), `y_length` (length of cells in y direction). Output is
flow width.
"""
function get_flow_width(ldd::UInt8, x_length::Real, y_length::Real)::Real
    # take into account non-square cells
    # if ldd is 8 or 2 use x_length
    # if ldd is 4 or 6 use y_length
    if ldd == 2 || ldd == 8
        x_length
    elseif ldd == 4 || ldd == 6
        y_length
    else
        (x_length * y_length) / hypot(x_length, y_length)
    end
end

"""
    get_surface_width(flow_width, flow_length, land_area, river_location)

Return the surface flow width. Input `flow_width` (flow width), `flow_length` (flow length),
`land_area` (area covered by land (excluding river coverage)) and `river_location` (river
cell, boolean). Output is surface flow width `surface_width`.
"""
function get_surface_width(
    flow_width::Real,
    flow_length::Real,
    land_area::Real,
    river_location::Bool,
)::Real
    surface_width = river_location ? land_area / flow_length : flow_width
    return surface_width
end

# 2.5x faster power method
"Faster method for exponentiation"
pow(x::Real, y::Real)::Real = exp(y * log(x))

function sum_at(A::AbstractVector{T}, inds::AbstractVector{Int})::T where {T}
    mapreduce(i -> A[i], +, inds; init = zero(T))
end

sum_at(f::Function, inds::AbstractVector{Int}; T::Type{<:Number} = Float64) =
    mapreduce(f, +, inds; init = zero(T))

# https://juliaarrays.github.io/StaticArrays.jl/latest/pages/api/#Arrays-of-static-arrays-1
function svectorscopy(x::Matrix{T}, ::Val{N})::Vector{SVector{N, T}} where {T, N}
    size(x, 1) == N || error("sizes mismatch")
    isbitstype(T) || error("use for bitstypes only")
    return copy(reinterpret(SVector{N, T}, vec(x)))
end

"""
    get_flow_fraction_to_river(graph, ldd, inds_river, slope)

Return flow `fraction` to a river cell (at index `j`) based on the ratio of the land surface
`slope` at index `j` to the sum of the land surface `slope` at index `j` and at river cell
index `i`.
"""
function get_flow_fraction_to_river(
    graph::SimpleDiGraph{Int},
    ldd::Vector{UInt8},
    inds_river::Vector{Int},
    slope::Vector{<:Real},
)::Vector{Float64}
    n = length(slope)
    fraction = zeros(n)
    for i in inds_river
        nbs = inneighbors(graph, i)
        for j in nbs
            if ldd[j] != ldd[i]
                fraction[j] = slope[j] / (slope[i] + slope[j])
            end
        end
    end
    return fraction
end

"""
    equal_size_vectors(x)

Used in the structs of arrays to ensure all vectors are of equal length.

`equal_size_vectors(([1,2], [1,2,3]))` would throw an ArgumentError.
`equal_size_vectors(([4,5], [4,5]))` would pass.
`equal_size_vectors((1, [4,5], [4,5]))` would also pass, since `1` is not an AbstractVector.
"""
function equal_size_vectors(x::Tuple)
    # all vectors in this struct should be the same size
    inds_vec = findall(arg -> isa(arg, AbstractVector), x)
    n = length(x[inds_vec[1]])
    x_vec = x[inds_vec]

    for arr in x_vec
        if length(arr) != n
            throw(ArgumentError("Not all vectors are of equal length"))
        end
    end
    return x
end

"""
    tosecond(x::Period)

Convert a Period into a Float64, which represents the number of seconds. Will fail if this
is not well defined, such as for Month.

# Examples
```julia-repl
julia> tosecond(Day(1))
86400.0
```
"""
tosecond(x::Hour) = Float64(Dates.value(Second(x)))
tosecond(x::Minute) = Float64(Dates.value(Second(x)))
tosecond(x::T) where {T <: DatePeriod} = Float64(Dates.value(Second(x)))
tosecond(x::T) where {T <: TimePeriod} = x / convert(T, Second(1))

"Return julian day of year (leap days are not counted)"
function julian_day(time::TimeType)::Int
    # for all years February 28 is day 59 and March 1 is day 60.
    day = dayofyear(time) - (isleapyear(time) && dayofyear(time) > 60)
    return day
end

"Partition indices with at least size `basesize`"
function _partition(xs::Integer, basesize::Integer)
    n = Int(max(1, xs ÷ basesize))
    return (Int(1 + ((i - 1) * xs) ÷ n):Int((i * xs) ÷ n) for i in 1:n)
end

"""
    threaded_foreach(f, x::AbstractArray; basesize::Integer)

Run function `f` in parallel by spawning tasks (nthreads <= 8), each task iterates over a
chunk of size `basesize`. For nthreads > 8 run function `f` in parallel with
`Polyester@batch` with `minbatch` equal to `basesize`.
"""
function threaded_foreach(f::Function, x::AbstractArray; basesize::Integer)::Nothing
    if Threads.nthreads() <= 8
        len = length(x)
        partitions = _partition(len, basesize)
        if length(partitions) > 1 && Threads.nthreads() > 1
            @sync for p in partitions
                Threads.@spawn begin
                    for i in eachindex(p)
                        f(@inbounds p[i])
                    end
                end
            end
        else
            for i in eachindex(x)
                f(@inbounds x[i])
            end
        end
    else
        @batch per = thread minbatch = basesize for i in eachindex(x)
            f(@inbounds x[i])
        end
    end
    return nothing
end

"""
    hydraulic_conductivity_at_depth(p::KvExponential, vertical_hydraulic_conductivity_factor, z, i, n)
    hydraulic_conductivity_at_depth(p::KvExponentialConstant, vertical_hydraulic_conductivity_factor, z, i, n)
    hydraulic_conductivity_at_depth(p::KvLayered, vertical_hydraulic_conductivity_factor, z, i, n)
    hydraulic_conductivity_at_depth(p::KvLayeredExponential, vertical_hydraulic_conductivity_factor, z, i, n)

Return vertical hydraulic conductivity `kv_z` at depth `z` for index `i` using multiplication
factor `kv_frac` at soil layer `n` and vertical hydraulic conductivity profile `p`.
"""
function hydraulic_conductivity_at_depth(
    p::KvExponential,
    vertical_hydraulic_conductivity_factor,
    z,
    i,
    n,
)
    kv_z =
        vertical_hydraulic_conductivity_factor[i][n] *
        p.kv_0[i] *
        exp(-p.hydraulic_conductivity_scale_parameter[i] * z)
    return kv_z
end

function hydraulic_conductivity_at_depth(
    p::KvExponentialConstant,
    vertical_hydraulic_conductivity_factor,
    z,
    i,
    n,
)
    (; kv_0, hydraulic_conductivity_scale_parameter) = p.exponential
    if z < p.z_exp[i]
        kv_z =
            vertical_hydraulic_conductivity_factor[i][n] *
            kv_0[i] *
            exp(-hydraulic_conductivity_scale_parameter[i] * z)
    else
        kv_z =
            vertical_hydraulic_conductivity_factor[i][n] *
            kv_0[i] *
            exp(-hydraulic_conductivity_scale_parameter[i] * p.z_exp[i])
    end
    return kv_z
end

function hydraulic_conductivity_at_depth(
    p::KvLayered,
    vertical_hydraulic_conductivity_factor,
    z,
    i,
    n,
)
    kv_z = vertical_hydraulic_conductivity_factor[i][n] * p.kv[i][n]
    return kv_z
end

function hydraulic_conductivity_at_depth(
    p::KvLayeredExponential,
    vertical_hydraulic_conductivity_factor,
    z,
    i,
    n,
)
    return if z < p.z_layered[i]
        vertical_hydraulic_conductivity_factor[i][n] * p.kv[i][n]
    else
        n = p.nlayers_kv[i]
        vertical_hydraulic_conductivity_factor[i][n] *
        p.kv[i][n] *
        exp(-p.hydraulic_conductivity_scale_parameter[i] * (z - p.z_layered[i]))
    end
end

"""
    bounded_divide(x, y; max = 1.0, default = 0.0)

Return the division of `x` by `y`, bounded by a maximum value `max`, when `y` > 0.0.
Otherwise return a `default` value.
"""
function bounded_divide(x::Real, y::Real; max::Real = 1.0, default::Real = 0.0)::Real
    z = y > 0.0 ? min(x / y, max) : default
    return z
end

"""
    bounded_power(base, power)

Computes min(base^power, 1) without computing the power
if the result is known to be larger than 1.
Assumes base, power > 0
"""
function bounded_power(base::T, power) where {T}
    return if base > 1
        one(T)
    else
        pow(base, power)
    end
end

"""
The sine of the slope in radians;
sin(arctan(x)) = x / √(1 + x²)
"""
sin_slope(slope) = slope / sqrt(1 + slope^2)

"Set lower bound for drainable porosity"
function lower_bound_drainable_porosity(theta_s, theta_fc; lower_bound = 0.02)
    return max(theta_s - theta_fc, lower_bound)
end
