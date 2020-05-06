const mv = NaN

Base.@kwdef struct LateralSSF{T,N}
    kh₀::Vector{T}                          # Horizontal hydraulic conductivity at soil surface [mm Δt⁻¹]
    f::Vector{T}                            # A scaling parameter [mm⁻¹] (controls exponential decline of kh₀)
    soilthickness::Vector{T}                # Soil thickness [mm]
    θₑ::Vector{T}                           # Effective porosity [-]
    Δt::Dates.Second                        # model time step [s]
    βₗ::Vector{T}                           # Slope [m m⁻¹]
    dl::Vector{T}                           # Drain length [mm]
    dw::Vector{T}                           # Flow width [mm]
    zi::Vector{T} = fill(mv, N)             # Pseudo-water table depth [mm] (top of the saturated zone)
    exfiltwater::Vector{T} = fill(mv, N)    # Exfiltration [mm]  (groundwater above surface level, saturated excess conditions)
    recharge::Vector{T} = fill(mv, N)       # Net recharge to saturated store [mm]
    ssf::Vector{T} =
        ((kh₀ .* βₗ) ./ f) .* (exp.(-f .* zi) - exp.(-f .* soilthickness)) .* dw    # Subsurface flow [mm³ Δt⁻¹]
    ssfmax::Vector{T} = ((kh₀ .* βₗ) ./ f) .* (1.0 .- exp.(-f .* soilthickness))     # Maximum subsurface flow [mm² Δt⁻¹]
end

"""
    statenames(::Type{LateralSSF})

Returns Array{Symbol,1} for extracting model state fields from SBM struct.
"""
function statenames(::Type{LateralSSF})

    # depends on ini file settings (optional: glaciers, snow, irrigation)
    states = [:ssf]
    #TODO: (warm) states read from netcdf file or cold state (reinit=1, setting in ini file)

end

function update(ssf::LateralSSF, dag, toposort, n)

    ssfn = zeros(n)
    zi = zeros(n)
    exfiltwater = zeros(n)

    for v in toposort
        upstream_nodes = inneighbors(dag, v)
        ssfin =
            isempty(upstream_nodes) ? 0.0 : sum(ssfn[i] for i in upstream_nodes)
        ssfn[v], zi[v], exfiltwater[v] = kinematic_wave_ssf(
            ssfin,
            ssf.ssf[v],
            ssf.zi[v],
            ssf.recharge[v],
            ssf.kh₀[v],
            ssf.βₗ[v],
            ssf.θₑ[v],
            ssf.f[v],
            ssf.soilthickness[v],
            Float64(ssf.Δt.value),
            ssf.dl[v],
            ssf.dw[v],
            ssf.ssfmax[v],
        )
    end
    return LateralSSF{Float64,n}(
        kh₀ = ssf.kh₀,
        f = ssf.f,
        soilthickness = ssf.soilthickness,
        θₑ = ssf.θₑ,
        Δt = ssf.Δt,
        βₗ = ssf.βₗ,
        dl = ssf.dl,
        dw = ssf.dw,
        zi = ssf.zi,
        exfiltwater = exfiltwater,
        recharge = ssf.recharge,
        ssf = ssfn,
    )

end
