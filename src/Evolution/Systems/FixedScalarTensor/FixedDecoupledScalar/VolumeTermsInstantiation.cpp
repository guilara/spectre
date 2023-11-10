// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <optional>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/VolumeTermsImpl.tpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/System.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/TimeDerivativeTerms.hpp"

namespace evolution::dg::Actions::detail {
template void volume_terms<::fe::DecoupledScalar::TimeDerivativeTerms>(
    const gsl::not_null<Variables<db::wrap_tags_in<
        ::Tags::dt,
        typename ::fe::DecoupledScalar::System::variables_tag::tags_list>>*>
        dt_vars_ptr,
    const gsl::not_null<Variables<db::wrap_tags_in<
        ::Tags::Flux, typename ::fe::DecoupledScalar::System::flux_variables,
        tmpl::size_t<3>, Frame::Inertial>>*>
        volume_fluxes,
    const gsl::not_null<Variables<db::wrap_tags_in<
        ::Tags::deriv,
        typename ::fe::DecoupledScalar::System::gradient_variables,
        tmpl::size_t<3>, Frame::Inertial>>*>
        partial_derivs,
    const gsl::not_null<
        Variables<typename ::fe::DecoupledScalar::System::
                      compute_volume_time_derivative_terms::temporary_tags>*>
        temporaries,
    const gsl::not_null<Variables<db::wrap_tags_in<
        ::Tags::div,
        db::wrap_tags_in<::Tags::Flux,
                         typename ::fe::DecoupledScalar::System::flux_variables,
                         tmpl::size_t<3>, Frame::Inertial>>>*>
        div_fluxes,
    const Variables<
        typename ::fe::DecoupledScalar::System::variables_tag::tags_list>&
        evolved_vars,
    const ::dg::Formulation dg_formulation, const Mesh<3>& mesh,
    [[maybe_unused]] const tnsr::I<DataVector, 3, Frame::Inertial>&
        inertial_coordinates,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>&
        logical_to_inertial_inverse_jacobian,
    [[maybe_unused]] const Scalar<DataVector>* const det_inverse_jacobian,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>& mesh_velocity,
    const std::optional<Scalar<DataVector>>& div_mesh_velocity,
    // GH argument variables
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::aa<DataVector, 3>& pi, const tnsr::iaa<DataVector, 3>& phi,
    const Scalar<DataVector>& gamma0, const Scalar<DataVector>& gamma1,
    const Scalar<DataVector>& gamma2,
    const ::gh::gauges::GaugeCondition& gauge_condition,
    const Mesh<3>& mesh_for_rhs, const double& time,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inverse_jacobian,
    const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
        mesh_velocity_gh,
    // Scalar argument variables
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3_st>& phi_scalar,
    // These appear with the same name as temporals for the other system
    const Scalar<DataVector>& lapse_scalar,
    const tnsr::I<DataVector, 3_st>& shift_scalar,

    const tnsr::i<DataVector, 3_st>& deriv_lapse,
    const tnsr::iJ<DataVector, 3_st>& deriv_shift,
    const tnsr::II<DataVector, 3_st>& upper_spatial_metric,
    const tnsr::I<DataVector, 3_st>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1_scalar,
    const Scalar<DataVector>& gamma2_scalar,
    // Scalar Tensor Extra Argument tags (Sources)
    const Scalar<DataVector>& scalar_source,

    // Scalar driver argument tags
    const Scalar<DataVector>& psi_scalar_driver,
    const Scalar<DataVector>& pi_scalar_driver,
    const tnsr::i<DataVector, 3_st>& phi_scalar_driver,

    const Scalar<DataVector>& gamma1_scalar_driver,
    const Scalar<DataVector>& gamma2_scalar_driver,
    const Scalar<DataVector>& scalar_driver_source,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter);
}  // namespace evolution::dg::Actions::detail
