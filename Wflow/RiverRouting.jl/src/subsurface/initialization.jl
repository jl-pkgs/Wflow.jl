"""
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, parameters::AbstractLandParameters, kh_profile::KhExponential)
    initialize_lateral_ssf_model!(subsurface_flow_model::LateralSSFModel, parameters::AbstractLandParameters, kh_profile::KhExponentialConstant)

Initialize lateral subsurface variables `q` and `q_max` using horizontal hydraulic
conductivity profile `kh_profile`.
"""
function initialize_lateral_ssf_model!(
    subsurface_flow_model::LateralSSFModel,
    parameters::AbstractLandParameters,
    kh_profile::KhExponential,
)
    (; kh_0, hydraulic_conductivity_scale_parameter) = kh_profile
    (; q, q_max, water_table_depth) = subsurface_flow_model.variables
    (; soil_thickness) = subsurface_flow_model.parameters
    (; slope, flow_width) = parameters

    @. q_max =
        ((kh_0 * slope) / hydraulic_conductivity_scale_parameter) *
        (1.0 - exp(-hydraulic_conductivity_scale_parameter * soil_thickness))
    @. q =
        ((kh_0 * slope) / hydraulic_conductivity_scale_parameter) *
        (
            exp(-hydraulic_conductivity_scale_parameter * water_table_depth) -
            exp(-hydraulic_conductivity_scale_parameter * soil_thickness)
        ) *
        flow_width
    return nothing
end

function initialize_lateral_ssf_model!(
    subsurface_flow_model::LateralSSFModel,
    parameters::AbstractLandParameters,
    kh_profile::KhExponentialConstant,
)
    (; kh_0, hydraulic_conductivity_scale_parameter) = kh_profile.exponential
    (; z_exp) = kh_profile
    (; q, q_max, water_table_depth) = subsurface_flow_model.variables
    (; soil_thickness) = subsurface_flow_model.parameters
    (; slope, flow_width) = parameters

    q_constant = @. kh_0 *
       exp(-hydraulic_conductivity_scale_parameter * z_exp) *
       slope *
       (soil_thickness - z_exp)
    for i in eachindex(q)
        q_max[i] =
            ((kh_0[i] * slope[i]) / hydraulic_conductivity_scale_parameter[i]) *
            (1.0 - exp(-hydraulic_conductivity_scale_parameter[i] * z_exp[i])) +
            q_constant[i]
        if water_table_depth[i] < z_exp[i]
            q[i] =
                (
                    ((kh_0[i] * slope[i]) / hydraulic_conductivity_scale_parameter[i]) * (
                        exp(
                            -hydraulic_conductivity_scale_parameter[i] *
                            water_table_depth[i],
                        ) - exp(-hydraulic_conductivity_scale_parameter[i] * z_exp[i])
                    ) + q_constant[i]
                ) * flow_width[i]
        else
            q[i] =
                kh_0[i] *
                exp(-hydraulic_conductivity_scale_parameter[i] * water_table_depth[i]) *
                slope[i] *
                (soil_thickness[i] - water_table_depth[i]) *
                flow_width[i]
        end
    end
    return nothing
end
