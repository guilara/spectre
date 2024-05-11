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
    const tnsr::I<DataVector, 3_st>& shift,
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
    const tnsr::aa<DataVector, 3>& tensor_driver);

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
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const Scalar<DataType>&,
                                    const Scalar<DataType>&,
                                    const Scalar<DataType>&,
                                    const Scalar<DataType>&) =
      &fe::ScalarDriver::Sources::compute_scalar_driver_source;
  using base = ScalarDriverSource;
};

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <typename Frame, typename DataType>
struct TargetScalarCompute : TargetScalar, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::WeylElectricScalar<DataType>,
                 gr::Tags::WeylMagneticScalar<DataType>,
                 CurvedScalarWave::Tags::Psi,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter,
                 ScalarTensor::Tags::ScalarMass>;
  using return_type = Scalar<DataType>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataType>&,
      const Scalar<DataType>&, const Scalar<DataType>&, const double,
      const double,
      const double) = &ScalarTensor::compute_scalar_curvature_source;
  using base = TargetScalar;
};

/*!
 * \brief Compute tag for the tensor driver source.
 */
struct TensorDriverSourceCompute
    : TensorDriverSource<DataVector, 3, Frame::Inertial>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, 3,
                                                            Frame::Inertial>,
                 fe::ScalarTensorDriver::Tags::TargetTensor<DataVector, 3,
                                                            Frame::Inertial>,
                 fe::ScalarTensorDriver::Tags::TauParameter,
                 fe::ScalarTensorDriver::Tags::SigmaParameter>;
  using return_type = tnsr::aa<DataVector, 3, Frame::Inertial>;
  static constexpr void (*function)(
      gsl::not_null<tnsr::aa<DataVector, 3>*> tensor_driver_source,
      const tnsr::aa<DataVector, 3>& tensor_driver,
      const tnsr::aa<DataVector, 3>& target_tensor,
      const Scalar<DataVector>& scalar_tau_parameter,
      const Scalar<DataVector>& scalar_sigma_parameter) =
      &fe::ScalarTensorDriver::Sources::compute_tensor_driver_source;
  using base = TensorDriverSource<DataVector, 3, Frame::Inertial>;
};

/*!
 * \brief Compute tag for the tensor driver target.
 * \details TODO : Replace with backreaction trace-reversed H tensor
 */
template <typename Frame, typename DataType>
struct TargetTensorCompute : TargetTensor<DataVector, 3, Frame::Inertial>,
                             db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, 3,
                                                            Frame::Inertial>>;
  using return_type = tnsr::aa<DataVector, 3, Frame::Inertial>;
  static constexpr void (*function)(
      gsl::not_null<tnsr::aa<DataVector, a>*> target_tensor,
      const tnsr::aa<DataVector, 3>& tensor_driver) =
      &fe::ScalarTensorDriver::Sources::compute_target_tensor;
  using base = TargetTensor<DataVector, 3, Frame::Inertial>;
};

}  // namespace fe::ScalarTensorDriver::Tags
