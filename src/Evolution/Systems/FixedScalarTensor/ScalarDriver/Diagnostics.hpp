// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "Options/String.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::ScalarDriver {

void driver_tracking_diagnostic(
    const gsl::not_null<Scalar<DataVector>*> diagnostic,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const Scalar<DataType>& scalar_tau_parameter,
    const Scalar<DataType>& scalar_sigma_parameter);

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
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>&,
                                    const Scalar<DataVector>&,
                                    const Scalar<DataVector>&,
                                    const Scalar<DataVector>&) =
      &fe::ScalarDriver::driver_tracking_diagnostic;
  using base = TrackingDiagnostic;
};

}  // namespace Tags

}  // namespace fe::ScalarDriver
