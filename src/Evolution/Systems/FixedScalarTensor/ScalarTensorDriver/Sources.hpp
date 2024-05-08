// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
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
    const tnsr::aa<DataVector>& tensor_driver,
    const tnsr::aa<DataVector, 3>& target_tensor,
    const Scalar<DataVector>& scalar_tau_parameter,
    const Scalar<DataVector>& scalar_sigma_parameter);

void compute_target_tensor(
    gsl::not_null<tnsr::aa<DataVector, a>*> target_tensor,
    const tnsr::aa<DataVector, 3>& tensor_driver);

}  // namespace fe::ScalarTensorDriver::Sources

namespace fe::ScalarTensorDriver::Tags {

/*!
 * \brief Compute tag for the tensor driver source.
 */
struct TensorDriverSourceCompute {};

/*!
 * \brief Compute tag for the tensor driver target.
 */
template <typename Frame, typename DataType>
struct TargetTensorCompute {};

}  // namespace fe::ScalarTensorDriver::Tags
