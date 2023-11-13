// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "Options/String.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::ScalarDriver {

void driver_tracking_diagnostic(
    const gsl::not_null<Scalar<DataVector>*> diagnostic,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter);

void shift_minus_mesh_velocity(
    const gsl::not_null<tnsr::I<DataVector, 3>*> result,
    const tnsr::I<DataVector, 3>& shift,
    const std::optional<tnsr::I<DataVector, 3>>& mesh_velocity);

namespace Tags {

/*!
 * \brief Compute a diagnostic for the scalar drive tracking of its target.
 *
 * \details Compute the difference between the scalar driver and its target.
 */
template <typename Frame, typename DataType>
struct TrackingDiagnosticCompute : TrackingDiagnostic, db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarDriver::Tags::Psi, fe::ScalarDriver::Tags::TargetPsi,
                 fe::ScalarDriver::Tags::ScalarTauParameter,
                 fe::ScalarDriver::Tags::ScalarSigmaParameter>;
  using return_type = Scalar<DataType>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataType>&,
      const Scalar<DataType>&, const Scalar<DataType>&,
      const Scalar<DataType>&) = &fe::ScalarDriver::driver_tracking_diagnostic;
  using base = TrackingDiagnostic;
};

/*!
 * \brief Compute the difference between the shift and mesh velocity.
 */
struct ShiftMinusMeshVelocityCompute : ShiftMinusMeshVelocity, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::Shift<DataVector, 3_st>,
                 Events::Tags::ObserverMeshVelocity<3_st, Frame::Inertial>>;
  using return_type = tnsr::I<DataVector, 3>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const tnsr::I<DataVector, 3>&,
      const std::optional<tnsr::I<DataVector, 3>>&) =
      &fe::ScalarDriver::shift_minus_mesh_velocity;
  using base = ShiftMinusMeshVelocity;
};

}  // namespace Tags

}  // namespace fe::ScalarDriver
