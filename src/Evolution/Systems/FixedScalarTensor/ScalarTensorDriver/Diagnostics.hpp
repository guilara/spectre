// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Diagnostics.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "Options/String.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::ScalarTensorDriver {

void tensor_driver_tracking_diagnostic(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> diagnostic,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::aa<DataVector, 3>& target_tensor,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter);

namespace Tags {

/*!
 * \brief Compute a diagnostic for the scalar driver tracking of its target.
 *
 * \details Compute the difference between the scalar driver and its target.
 */
template <typename Frame, typename DataType>
struct ScalarTrackingDiagnosticCompute : ScalarTrackingDiagnostic,
                                         db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarTensorDriver::Tags::Psi,
                 fe::ScalarTensorDriver::Tags::TargetScalar,
                 fe::ScalarTensorDriver::Tags::TauParameter,
                 fe::ScalarTensorDriver::Tags::SigmaParameter>;
  using return_type = Scalar<DataType>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataType>*> diagnostic,
      const Scalar<DataType>& psi, const Scalar<DataType>& target_psi,
      const Scalar<DataType>& scalar_tau_parameter,
      const Scalar<DataType>& scalar_sigma_parameter) =
      &fe::ScalarDriver::driver_tracking_diagnostic;
  using base = ScalarTrackingDiagnostic;
};

/*!
 * \brief Compute a diagnostic for the tensor driver tracking of its target.
 *
 * \details Compute the difference between the tensor driver and its target.
 */
template <typename Frame, typename DataType>
struct TensorTrackingDiagnosticCompute
    : TensorTrackingDiagnostic<DataType, 3, Frame>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarTensorDriver::Tags::TensorDriver<DataType, 3, Frame>,
                 fe::ScalarTensorDriver::Tags::TargetTensor<DataType, 3, Frame>,
                 fe::ScalarTensorDriver::Tags::TauParameter,
                 fe::ScalarTensorDriver::Tags::SigmaParameter>;
  using return_type = tnsr::aa<DataType, 3, Frame>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataType, 3, Frame>*> diagnostic,
      const tnsr::aa<DataType, 3, Frame>& tensor_driver,
      const tnsr::aa<DataType, 3, Frame>& target_tensor,
      const Scalar<DataType>& scalar_tau_parameter,
      const Scalar<DataType>& scalar_sigma_parameter) =
      &tensor_driver_tracking_diagnostic;
  using base = TensorTrackingDiagnostic<DataType, 3, Frame>;
};

}  // namespace Tags

}  // namespace fe::ScalarTensorDriver
