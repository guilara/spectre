// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <optional>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/VolumeTermsImpl.tpp"
#include "Evolution/Systems/ScalarTensor/TwoScalars/System.hpp"
#include "Evolution/Systems/ScalarTensor/TwoScalars/TimeDerivative.hpp"

namespace evolution::dg::Actions::detail {
template void volume_terms<::ScalarTensor::TwoScalars::TimeDerivative>(
    const gsl::not_null<Variables<
        db::wrap_tags_in<::Tags::dt, typename ::ScalarTensor::TwoScalars::
                                         System::variables_tag::tags_list>>*>
        dt_vars_ptr,
    const gsl::not_null<Variables<db::wrap_tags_in<
        ::Tags::Flux,
        typename ::ScalarTensor::TwoScalars::System::flux_variables,
        tmpl::size_t<3>, Frame::Inertial>>*>
        volume_fluxes,
    const gsl::not_null<Variables<db::wrap_tags_in<
        ::Tags::deriv,
        typename ::ScalarTensor::TwoScalars::System::gradient_variables,
        tmpl::size_t<3>, Frame::Inertial>>*>
        partial_derivs,
    const gsl::not_null<
        Variables<typename ::ScalarTensor::TwoScalars::System::
                      compute_volume_time_derivative_terms::temporary_tags>*>
        temporaries,
    const gsl::not_null<Variables<db::wrap_tags_in<
        ::Tags::div, db::wrap_tags_in<::Tags::Flux,
                                      typename ::ScalarTensor::TwoScalars::
                                          System::flux_variables,
                                      tmpl::size_t<3>, Frame::Inertial>>>*>
        div_fluxes,
    const Variables<
        typename ::ScalarTensor::TwoScalars::System::variables_tag::tags_list>&
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
    const tnsr::i<DataVector, 3>& phi_scalar,
    const Scalar<DataVector>& lapse_scalar,
    const tnsr::I<DataVector, 3>& shift_scalar,
    const tnsr::i<DataVector, 3>& deriv_lapse,
    const tnsr::iJ<DataVector, 3>& deriv_shift,
    const tnsr::II<DataVector, 3>& upper_spatial_metric,
    const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1_scalar,
    const Scalar<DataVector>& gamma2_scalar,
    // Scalar argument variables
    const Scalar<DataVector>& pi_scalar_2,
    const tnsr::i<DataVector, 3>& phi_scalar_2,
    const Scalar<DataVector>& gamma1_scalar_2,
    const Scalar<DataVector>& gamma2_scalar_2,

    // Scalar Tensor extra argument tags
    const Scalar<DataVector>& scalar_source,
    const Scalar<DataVector>& scalar_source_2);
}  // namespace evolution::dg::Actions::detail
