const MISSING_VALUE = Float64(NaN)
const GRAVITATIONAL_ACCELERATION = 9.80665 # m s⁻²
const SH = NamedTuple{(:H, :S), Tuple{Vector{Float64}, Vector{Float64}}}
const HQ = NamedTuple{(:H, :Q), Tuple{Vector{Float64}, Matrix{Float64}}}

"Faster method for exponentiation."
pow(x::Real, y::Real)::Real = exp(y * log(x))

"Sigmoid S-shaped curve."
scurve(x::Real, a::Real, b::Real, c::Real)::Real = inv(b + exp(-c * (x - a)))

function sum_at(values::AbstractVector{T}, indices::AbstractVector{Int})::T where {T}
    return mapreduce(index -> values[index], +, indices; init = zero(T))
end

sum_at(func::Function, indices::AbstractVector{Int}; T::Type{<:Number} = Float64) =
    mapreduce(func, +, indices; init = zero(T))

function to_enumx(enum_type, value::Int)
    options = instances(enum_type)
    if 1 ≤ value ≤ length(options)
        return options[value]
    end
    options_repr = repr(MIME("text/plain"), enum_type)
    error(
        "Cannot convert $value to $enum_type, there are only $(length(options)) options:\n$options_repr.",
    )
end

"""
    set_layerthickness(reference_depth, cumulative_depth, thickness)

Calculate active soil-layer thicknesses down to `reference_depth`. Depths and thicknesses
are in metres.
"""
function set_layerthickness(
    reference_depth::Real,
    cumulative_depth::SVector,
    thickness::SVector{N, Float64},
)::SVector{N, Float64} where {N}
    active_thickness = thickness .* MISSING_VALUE
    for layer_idx in eachindex(active_thickness)
        if reference_depth > cumulative_depth[layer_idx + 1]
            active_thickness = setindex(active_thickness, thickness[layer_idx], layer_idx)
        elseif reference_depth - cumulative_depth[layer_idx] > 0.0
            active_thickness = setindex(
                active_thickness,
                reference_depth - cumulative_depth[layer_idx],
                layer_idx,
            )
        end
    end
    return active_thickness
end

"Return the number of non-NaN active soil layers."
number_of_active_layers(thickness::SVector)::Int = length(thickness) - sum(isnan, thickness)

"Set a lower bound for drainable porosity."
lower_bound_drainable_porosity(theta_s, theta_fc; lower_bound = 0.02) =
    max(theta_s - theta_fc, lower_bound)

"Return Julian day of year while omitting leap days."
function julian_day(time::TimeType)::Int
    ordinal_day = dayofyear(time)
    return ordinal_day - (isleapyear(time) && ordinal_day > 60)
end

"Read a storage-water-level curve from CSV."
function read_sh_csv(path)
    data, header = readdlm(path, ',', Float64; header = true)
    names = Symbol.(vec(header))
    return NamedTuple{Tuple(names)}(
        Tuple(data[:, column_idx] for column_idx in axes(data, 2)),
    )
end

"Read a water-level-discharge curve from CSV."
function read_hq_csv(path)
    data = readdlm(path, ',', Float64; skipstart = 1)
    return (; H = data[:, 1], Q = data[:, 2:end])
end

"Partition indices with at least size `basesize`."
function _partition(length_x::Integer, basesize::Integer)
    n_partitions = Int(max(1, length_x ÷ basesize))
    return (
        Int(1 + ((partition_idx - 1) * length_x) ÷ n_partitions):Int(
            (partition_idx * length_x) ÷ n_partitions,
        ) for partition_idx in 1:n_partitions
    )
end

"Run `func` over an array using Julia tasks or Polyester threads."
function threaded_foreach(func::Function, values::AbstractArray; basesize::Integer)::Nothing
    if Threads.nthreads() <= 8
        partitions = _partition(length(values), basesize)
        if length(partitions) > 1 && Threads.nthreads() > 1
            @sync for partition in partitions
                Threads.@spawn for partition_idx in eachindex(partition)
                    func(@inbounds partition[partition_idx])
                end
            end
        else
            for value_idx in eachindex(values)
                func(@inbounds values[value_idx])
            end
        end
    else
        @batch per = thread minbatch = basesize for value_idx in eachindex(values)
            func(@inbounds values[value_idx])
        end
    end
    return nothing
end
