// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/WeylTensor.hpp"

#include <cstddef>

#include "DataStructures/LeviCivitaIterator.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/VectorImpl.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"

namespace gr {

template <size_t SpatialDim, typename Frame, typename DataType>
void weyl_tensor_from_electric_magnetic(
    gsl::not_null<tnsr::abcc<DataType, VolumeDim, Frame>*> weyl_result){
  destructive_resize_components(weyl_result,
                                get_size(get<0, 0>(inverse_spatial_metric)));
  *weyl_result = make_with_value<tnsr::abcc<DataType, VolumeDim, Frame>>(
      get<0, 0>(inverse_spatial_metric), 0.0);

  get(*weyl_result) +=
      4 * weyl_electric.get() *
      (spatial_metric.get() + normal_vector.get() * normal_vector.get());
  get(*weyl_result) +=
      -1.0 * /* 3 levi-civita */ normal_vector.get() * weyl_magnetic.get();

  // Then apply projector to impose symmetries
}

template <size_t SpatialDim, typename Frame, typename DataType>
tnsr::abcc<DataType, VolumeDim, Frame> weyl_tensor_from_electric_magnetic() {
//
}

template <typename DataType, typename Frame>
void weyl_square_scalar_compute(
    gsl::not_null<Scalar<DataType>*> weyl_square_result,
    tnsr::abcc<DataType, 4_st, Frame>& weyl_tensor,
    tnsr::AA<DataType, 4_st, Frame>& inverse_spacetime_metric) {
    //
    destructive_resize_components(weyl_square_result,
                                get_size(get<0, 0>(inverse_spacetime_metric)));
    *weyl_square_result = make_with_value<Scalar<DataType>>(
        get<0, 0>(inverse_spacetime_metric), 0.0);

    for (size_t a1 = 0; a1 < 4; ++a1) {
    for (size_t b1 = 0; b1 < 4; ++b1) {
        for (size_t c1 = 0; c1 < 4; ++c1) {
        for (size_t d1 = 0; d1 < 4; ++d1) {
            for (size_t a2 = 0; a2 < 4; ++a2) {
            for (size_t b2 = 0; b2 < 4; ++b2) {
                for (size_t c2 = 0; c2 < 4; ++c2) {
                for (size_t d2 = 0; d2 < 4; ++d2) {
                    get(*weyl_square_result) +=
                        weyl_tensor.get(a1, b1, c1, d1) *
                        inverse_spacetime_metric.get(a1, a2) *
                        inverse_spacetime_metric.get(b1, b2) *
                        inverse_spacetime_metric.get(c1, c2) *
                        inverse_spacetime_metric.get(d1, d2) *
                        weyl_tensor.get(a2, b2, c2, d2);
                }
                }
            }
            }
        }
        }
    }
    }
}

template <typename DataType, typename Frame>
Scalar<DataType> weyl_square_scalar_compute(
    tnsr::abcc<DataType, 4_st, Frame>& weyl_tensor,
    tnsr::AA<DataType, 4_st, Frame>& inverse_spacetime_metric) {
Scalar<DataType> weyl_square_result{get<0, 0>(inverse_spacetime_metric)};
weyl_square_scalar_compute<DataType, Frame>(
    make_not_null(&weyl_square_result), weyl_tensor, inverse_spacetime_metric);
return weyl_square_result;
}

}  // namespace gr
