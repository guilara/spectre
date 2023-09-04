// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/ScalarTensor/GBTensorSourceTerms.hpp"

#include "DataStructures/LeviCivitaIterator.hpp"

namespace ScalarTensor {

void gb_H_tensor(const gsl::not_null<tnsr::aa<DataVector, 3, Frame>*> result,
                 const tnsr::ii<DataVector, 3, Frame>& weyl_electric,
                 const tnsr::ii<DataVector, 3, Frame>& weyl_magnetic,
                 const tnsr::ii<DataVector, 3>& spatial_metric,
                 const tnsr::II<DataVector, 3>& inverse_spatial_metric,
                 const Scalar<DataVector>& sqrt_det_spatial_metric const
                     tnsr::aa<DataVector, Dim>& spacetime_metric,
                 const tnsr::A<DataVector, 3>& normal_spacetime_vector,
                 const tnsr::a<DataVector, 3>& normal_spacetime_one_form,
                 const tnsr::aa<DataVector, 3>& dd_coupling_function) {
  // Raise one index of the magnetic part with the spatial metric
  weyl_magnetic_down_up.get(i, j) =
      weyl_magnetic.get(i, k) * inverse_spatial_metric.get(k, j);

  tnsr::abcd<DataVector, 3> spacetime_weyl_I =
      make_with_value<tnsr::abcd<DataVector, 3>>(get<0, 0>(spatial_metric),
                                                 0.0);
  // Electric part
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      for (size_t c = 0; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          if (a == 0 or b == 0) {
            spacetime_weyl_I.get(a, b, c, d) += 0.0;
          } else if (c > 0 and d > 0) {
            i = a - 1;
            j = b - 1;
            k = c - 1;
            l = d - 1;
            spacetime_weyl_I.get(a, b, c, d) +=
                weyl_electric.get(i, j) *
                (spatial_metric.get(k, l) +
                 normal_spacetime_one_form.get(k) *
                     normal_spacetime_one_form.get(l));
          } else {
            i = a - 1;
            j = b - 1;
            spacetime_weyl_I.get(a, b, c, d) +=
                weyl_electric.get(i, j) * (normal_spacetime_one_form.get(c) *
                                           normal_spacetime_one_form.get(d));
          }
        }
      }
    }
  }

  // Magnetic part
  for (size_t a = 1; a < 4; ++a) {
    for (size_t b = 1; b < 4; ++b) {
      for (size_t c = 1; c < 4; ++c) {
        for (size_t d = 0; d < 4; ++d) {
          i = a - 1;
          j = b - 1;
          k = c - 1;
          for (size_t m = 0; m < 3; ++m) {
            spacetime_weyl_I.get(a, b, c, d) +=
                three_epsilon.get(i, j, m) * normal_spacetime_one_form.get(d) *
                weyl_magnetic_down_up.get(k, m);
          }
        }
      }
    }
  }
}

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
    const tnsr::aa<DataVector, 3>& dd_coupling_function) {
  tnsr::abcd<DataVector, 3> spacetime_weyl_I =
      make_with_value<tnsr::abcd<DataVector, 3>>(get<0, 0>(spatial_metric),
                                                 0.0);
  // Electric part
}

}  // namespace ScalarTensor
