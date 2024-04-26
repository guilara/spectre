// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor {

// template <typename Frame>
void DDKG_tensor_from_projections(
    const gsl::not_null<tnsr::aa<DataVector, 3>*> DDKG_tensor_result,

    // Metric quantities
    const Scalar<DataVector>& lapse, const tnsr::I<DataVector, 3>& shift,
    const tnsr::aa<DataVector, 3>& spacetime_metric,
    const tnsr::II<DataVector, 3>& inverse_spatial_metric,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature,
    const tnsr::Ijj<DataVector, 3>& spatial_christoffel_second_kind,

    // Scalar quantities
    const Scalar<DataVector>& psi_scalar, const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,

    // Scalar gradients
    const tnsr::i<DataVector, 3>& d_psi_scalar,
    const tnsr::i<DataVector, 3>& d_pi_scalar,
    const tnsr::ij<DataVector, 3>& d_phi_scalar,

    // Provide them with RHS compute tags
    const Scalar<DataVector>& dt_psi_scalar,
    const Scalar<DataVector>& dt_pi_scalar,
    const tnsr::i<DataVector, 3>& dt_phi_scalar);

/*
void order_reduced_gb_H_tensor(
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
*/

}  // namespace ScalarTensor
