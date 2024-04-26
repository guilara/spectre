// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

namespace ScalarTensor {

void gb_H_tensor(const gsl::not_null<tnsr::aa<DataVector, 3, Frame>*> result,
                 const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
                 const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
                 const tnsr::ii<DataVector, 3>& spatial_metric,
                 const tnsr::II<DataVector, 3>& inverse_spatial_metric,
                 const Scalar<DataVector>& sqrt_det_spatial_metric
                 const tnsr::aa<DataVector, Dim>& spacetime_metric,
                 const tnsr::A<DataVector, 3>& normal_spacetime_vector,
                 const tnsr::a<DataVector, 3>& normal_spacetime_one_form,
                 const tnsr::aa<DataVector, 3>& dd_coupling_function);

void gb_H_tensor_with_tenex(
    const gsl::not_null<tnsr::aa<DataVector, 3, Frame>*> result,
    const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
    const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const Scalar<DataVector>& sqrt_det_spatial_metric const
        tnsr::aa<DataVector, Dim>& spacetime_metric,
    const tnsr::A<DataVector, 3>& normal_spacetime_vector,
    const tnsr::a<DataVector, 3>& normal_spacetime_one_form,
    const tnsr::aa<DataVector, 3>& dd_coupling_function);

}  // namespace ScalarTensor
