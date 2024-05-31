// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Sources.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::ScalarTensorDriver::Sources {

void add_tensor_driver_source_to_dt_pi(
    gsl::not_null<tnsr::aa<DataVector, 3>*> dt_pi,
    const tnsr::aa<DataVector, 3>& tensor_driver_source,
    const Scalar<DataVector>& lapse);

void add_tensor_driver_friction_term_to_dt_pi(
    gsl::not_null<tnsr::aa<DataVector, 3>*> dt_pi,
    const tnsr::aa<DataVector, 3>& pi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter);

void compute_tensor_driver_source(
    gsl::not_null<tnsr::aa<DataVector, 3>*> tensor_driver_source,
    const tnsr::aa<DataVector, 3>& tensor_driver,
    const tnsr::aa<DataVector, 3>& target_tensor,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter);

void compute_target_tensor(
    gsl::not_null<tnsr::aa<DataVector, 3>*> target_tensor,
    const tnsr::aa<DataVector, 3>& trace_reversed_stress_energy);

}  // namespace fe::ScalarTensorDriver::Sources

namespace fe::ScalarTensorDriver::Tags {

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <typename Frame, typename DataType>
struct ScalarDriverSourceCompute : ScalarDriverSource, db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarTensorDriver::Tags::Psi,
                 fe::ScalarTensorDriver::Tags::TargetScalar,
                 fe::ScalarTensorDriver::Tags::TauParameter,
                 fe::ScalarTensorDriver::Tags::SigmaParameter>;
  using return_type = Scalar<DataType>;
  static constexpr void (*function)(
      gsl::not_null<Scalar<DataType>*> result, const Scalar<DataType>& psi,
      const Scalar<DataType>& target_psi,
      const Scalar<DataType>& scalar_tau_parameter,
      const Scalar<DataType>& scalar_sigma_parameter) =
      &fe::ScalarDriver::Sources::compute_scalar_driver_source;
  using base = ScalarDriverSource;
};

// /*!
//  * \brief Compute tag for the scalar driver source.
//  *
//  * \details Call compute_scalar_driver_source.
//  */
// template <typename Frame, typename DataType>
// struct TargetScalarCompute : TargetScalar, db::ComputeTag {
//   using argument_tags =
//       tmpl::list<gr::Tags::WeylElectricScalar<DataType>,
//                  gr::Tags::WeylMagneticScalar<DataType>,
//                  CurvedScalarWave::Tags::Psi,
//                  ScalarTensor::Tags::ScalarFirstCouplingParameter,
//                  ScalarTensor::Tags::ScalarSecondCouplingParameter,
//                  ScalarTensor::Tags::ScalarMass>;
//   using return_type = Scalar<DataType>;
//   static constexpr void (*function)(
//       gsl::not_null<Scalar<DataType>*> scalar_source,
//       const Scalar<DataType>& weyl_electric_scalar,
//       const Scalar<DataType>& weyl_magnetic_scalar, const Scalar<DataType>&
//       psi, const double first_coupling_psi, const double second_coupling_psi,
//       const double mass_psi) =
//       &ScalarTensor::compute_scalar_curvature_source;
//   using base = TargetScalar;
// };

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <typename Frame, typename DataType>
struct TargetScalarCompute : TargetScalar, db::ComputeTag {
  using argument_tags =
      tmpl::list<ScalarTensor::Tags::OrderReducedGBScalar,
                 CurvedScalarWave::Tags::Psi,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter,
                 ScalarTensor::Tags::ScalarMass>;
  using return_type = Scalar<DataType>;
  static constexpr void (*function)(
      gsl::not_null<Scalar<DataType>*> scalar_source,
      const Scalar<DataType>& order_reduced_gb_scalar,
      const Scalar<DataType>& psi, const double first_coupling_psi,
      const double second_coupling_psi, const double mass_psi) =
      &ScalarTensor::compute_order_reduced_scalar_curvature_source;
  using base = TargetScalar;
};

/*!
 * \brief Compute tag for the tensor driver source.
 */
template <typename Frame, typename DataType>
struct TensorDriverSourceCompute : TensorDriverSource<DataType, 3, Frame>,
                                   db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarTensorDriver::Tags::TensorDriver<DataType, 3, Frame>,
                 fe::ScalarTensorDriver::Tags::TargetTensor<DataType, 3, Frame>,
                 fe::ScalarTensorDriver::Tags::TauParameter,
                 fe::ScalarTensorDriver::Tags::SigmaParameter>;
  using return_type = tnsr::aa<DataType, 3, Frame>;
  static constexpr void (*function)(
      gsl::not_null<tnsr::aa<DataType, 3, Frame>*> tensor_driver_source,
      const tnsr::aa<DataType, 3, Frame>& tensor_driver,
      const tnsr::aa<DataType, 3, Frame>& target_tensor,
      const Scalar<DataType>& scalar_tau_parameter,
      const Scalar<DataType>& scalar_sigma_parameter) =
      &fe::ScalarTensorDriver::Sources::compute_tensor_driver_source;
  using base = TensorDriverSource<DataType, 3, Frame>;
};

// /*!
//  * \brief Compute tag for the tensor driver target.
//  * \details TODO : Replace with backreaction trace-reversed H tensor
//  */
// template <typename Frame, typename DataType>
// struct TargetTensorCompute : TargetTensor<DataType, 3, Frame>, db::ComputeTag
// {
//   using argument_tags = tmpl::list<
//       ScalarTensor::Tags::TraceReversedStressEnergy<DataType, 3, Frame>>;
//   using return_type = tnsr::aa<DataType, 3, Frame>;
//   static constexpr void (*function)(
//       gsl::not_null<tnsr::aa<DataType, 3, Frame>*> target_tensor,
//       const tnsr::aa<DataType, 3, Frame>& tensor_driver) =
//       &fe::ScalarTensorDriver::Sources::compute_target_tensor;
//   using base = TargetTensor<DataType, 3, Frame>;
// };

/*!
 * \brief Compute tag for the tensor driver target.
 * \details TODO : Replace with backreaction trace-reversed H tensor
 */
template <typename Frame, typename DataType>
struct TargetTensorCompute : TargetTensor<DataType, 3, Frame>, db::ComputeTag {
  using argument_tags =
      tmpl::list<ScalarTensor::Tags::OrderReducedTraceReversedStressEnergy>;
  using return_type = tnsr::aa<DataType, 3, Frame>;
  static constexpr void (*function)(
      gsl::not_null<tnsr::aa<DataType, 3, Frame>*> target_tensor,
      const tnsr::aa<DataType, 3, Frame>& tensor_driver) =
      &fe::ScalarTensorDriver::Sources::compute_target_tensor;
  using base = TargetTensor<DataType, 3, Frame>;
};

}  // namespace fe::ScalarTensorDriver::Tags
