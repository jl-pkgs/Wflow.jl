"""
    kh_layered_profile!(soil_model::SbmSoilModel, subsurface_flow_model::LateralSSFModel, kv_profile::KvLayered, dt)
    kh_layered_profile!(soil_model::SbmSoilModel, subsurface_flow_model::LateralSSFModel, kv_profile::KvLayeredExponential, dt)

Compute equivalent horizontal hydraulic conductivity `kh` [m d⁻¹] using vertical hydraulic
conductivity profile `kv_profile`.
"""
function kh_layered_profile!(
    soil_model::SbmSoilModel,
    subsurface_flow_model::LateralSSFModel,
    kv_profile::KvLayered,
)
    (; number_of_layers, cumulative_layer_depth, actual_layer_thickness, soil_thickness) =
        soil_model.parameters
    (; n_unsatlayers, water_table_depth) = soil_model.variables
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; horizontal_to_vertical_hydraulic_conductivity_ratio) =
        subsurface_flow_model.parameters

    for i in eachindex(kh)
        m = number_of_layers[i]

        if soil_thickness[i] > water_table_depth[i]
            transmissivity = 0.0
            _sumlayers = @view cumulative_layer_depth[i][2:end]
            n = max(n_unsatlayers[i], 1)
            transmissivity += (_sumlayers[n] - water_table_depth[i]) * kv_profile.kv[i][n]
            n += 1
            while n <= m
                transmissivity += actual_layer_thickness[i][n] * kv_profile.kv[i][n]
                n += 1
            end
            kh[i] =
                (transmissivity / (soil_thickness[i] - water_table_depth[i])) *
                horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        else
            kh[i] =
                kv_profile.kv[i][m] * horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        end
    end
    return nothing
end

function kh_layered_profile!(
    soil_model::SbmSoilModel,
    subsurface_flow_model::LateralSSFModel,
    kv_profile::KvLayeredExponential,
)
    (; number_of_layers, cumulative_layer_depth, actual_layer_thickness, soil_thickness) =
        soil_model.parameters
    (; nlayers_kv, z_layered, kv, hydraulic_conductivity_scale_parameter) = kv_profile
    (; n_unsatlayers, water_table_depth) = soil_model.variables
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; horizontal_to_vertical_hydraulic_conductivity_ratio) =
        subsurface_flow_model.parameters

    for i in eachindex(kh)
        m = number_of_layers[i]

        if soil_thickness[i] > water_table_depth[i]
            transmissivity = 0.0
            n = max(n_unsatlayers[i], 1)
            if water_table_depth[i] >= z_layered[i]
                zt = soil_thickness[i] - z_layered[i]
                j = nlayers_kv[i]
                transmissivity +=
                    kv[i][j] / hydraulic_conductivity_scale_parameter[i] * (
                        exp(
                            -hydraulic_conductivity_scale_parameter[i] *
                            (water_table_depth[i] - z_layered[i]),
                        ) - exp(-hydraulic_conductivity_scale_parameter[i] * zt)
                    )
                n = m
            else
                _sumlayers = @view cumulative_layer_depth[i][2:end]
                transmissivity += (_sumlayers[n] - water_table_depth[i]) * kv[i][n]
            end
            n += 1
            while n <= m
                if n > nlayers_kv[i]
                    zt = soil_thickness[i] - z_layered[i]
                    j = nlayers_kv[i]
                    transmissivity +=
                        kv[i][j] / hydraulic_conductivity_scale_parameter[i] *
                        (1.0 - exp(-hydraulic_conductivity_scale_parameter[i] * zt))
                    n = m
                else
                    transmissivity += actual_layer_thickness[i][n] * kv[i][n]
                end
                n += 1
            end
            kh[i] =
                (transmissivity / (soil_thickness[i] - water_table_depth[i])) *
                horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        else
            if water_table_depth[i] >= z_layered[i]
                j = nlayers_kv[i]
                kh[i] =
                    kv[i][j] *
                    exp(
                        -hydraulic_conductivity_scale_parameter[i] *
                        (water_table_depth[i] - z_layered[i]),
                    ) *
                    horizontal_to_vertical_hydraulic_conductivity_ratio[i]
            else
                kh[i] = kv[i][m] * horizontal_to_vertical_hydraulic_conductivity_ratio[i]
            end
        end
    end
    return nothing
end

kh_layered_profile!(
    soil_model::SbmSoilModel,
    subsurface_flow_model::LateralSSFModel,
    kv_profile::Union{KvExponential, KvExponentialConstant},
) = nothing

"""
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, soil_model::SbmSoilModel, parameters::LandParameters, kv_profile::KvLayered, dt)
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, soil_model::SbmSoilModel, parameters::LandParameters, kv_profile::KvLayeredExponential, dt)

Initialize lateral subsurface variables `q` and `q_max` using  vertical hydraulic
conductivity profile `kv_profile`.
"""
function initialize_lateral_ssf_model!(
    subsurface_flow_model::LateralSSFModel,
    soil_model::SbmSoilModel,
    parameters::LandParameters,
    kv_profile::KvLayered,
    dt,
)
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; number_of_layers, actual_layer_thickness) = soil_model.parameters
    (; q, q_max, water_table_depth) = subsurface_flow_model.variables
    (; horizontal_to_vertical_hydraulic_conductivity_ratio, soil_thickness) =
        subsurface_flow_model.parameters
    (; slope, flow_width) = parameters

    kh_layered_profile!(soil_model, subsurface_flow_model, kv_profile)
    for i in eachindex(q)
        q[i] = kh[i] * (soil_thickness[i] - water_table_depth[i]) * slope[i] * flow_width[i]
        kh_max = 0.0
        for j in 1:number_of_layers[i]
            kh_max += kv_profile.kv[i][j] * actual_layer_thickness[i][j]
        end
        kh_max *= horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        q_max[i] = kh_max * slope[i]
    end
    return nothing
end

function initialize_lateral_ssf_model!(
    subsurface_flow_model::LateralSSFModel,
    soil_model::SbmSoilModel,
    parameters::LandParameters,
    kv_profile::KvLayeredExponential,
    dt,
)
    (; q, q_max, water_table_depth) = subsurface_flow_model.variables
    (; horizontal_to_vertical_hydraulic_conductivity_ratio, soil_thickness) =
        subsurface_flow_model.parameters
    (; slope, flow_width) = parameters
    (; number_of_layers, actual_layer_thickness) = soil_model.parameters
    (; kh) = subsurface_flow_model.parameters.kh_profile
    (; kv, hydraulic_conductivity_scale_parameter, nlayers_kv, z_layered) = kv_profile

    kh_layered_profile!(soil_model, subsurface_flow_model, kv_profile)
    for i in eachindex(q)
        q[i] = kh[i] * (soil_thickness[i] - water_table_depth[i]) * slope[i] * flow_width[i]
        kh_max = 0.0
        for j in 1:number_of_layers[i]
            if j <= nlayers_kv[i]
                kh_max += kv[i][j] * actual_layer_thickness[i][j]
            else
                zt = soil_model.parameters.soil_thickness[i] - z_layered[i]
                k = max(j - 1, 1)
                kh_max +=
                    kv[i][k] / hydraulic_conductivity_scale_parameter[i] *
                    (1.0 - exp(-hydraulic_conductivity_scale_parameter[i] * zt))
                break
            end
        end
        kh_max = kh_max * horizontal_to_vertical_hydraulic_conductivity_ratio[i]
        q_max[i] = kh_max * slope[i]
    end
    return nothing
end

"""
Return water table change `dh` and exfiltration rate `exfilt`. For a falling water table
`dh` is based on subsurface net flux `net_flux` and specific yield `specific_yield`. For a
rising water table `dh` is based on `net_flux` and the unsaturated store capacity (per soil
layer). For a rising water table a dynamic specific yield is computed.
"""
function water_table_change(
    soil_model::SbmSoilModel,
    net_flux::Float64,
    specific_yield::Float64,
    i::Int,
    dt::Float64,
)
    (; n_unsatlayers, unsaturated_layer_thickness, unsaturated_layer_depth) =
        soil_model.variables
    (; theta_s, theta_r) = soil_model.parameters

    # effective porosity (difference between saturated and residual water content)
    theta_e = theta_s[i] - theta_r[i]

    if net_flux <= 0.0
        dh = net_flux * dt / specific_yield
    else
        dh = 0.0
        for k in n_unsatlayers[i]:-1:1
            capacity =
                max(
                    unsaturated_layer_thickness[i][k] * theta_e -
                    unsaturated_layer_depth[i][k],
                    0.0,
                ) / dt
            flux_layer = min(net_flux, capacity)
            if capacity <= net_flux
                # if unsaturated layer is fully saturated dh equals layer thickness
                dh += unsaturated_layer_thickness[i][k]
            else
                sy =
                    theta_e -
                    (unsaturated_layer_depth[i][k] / unsaturated_layer_thickness[i][k])
                dh += flux_layer * dt / sy
            end
            net_flux -= flux_layer
            net_flux == 0.0 && break
        end
    end
    exfilt = max(net_flux, 0.0)
    return dh, exfilt
end
