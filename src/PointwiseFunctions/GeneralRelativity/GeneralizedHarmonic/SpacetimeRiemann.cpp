// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpacetimeRiemann.hpp"

#include "DataStructures/DataVector.hpp"  // IWYU pragma: keep
#include "DataStructures/Tags/TempTensor.hpp"
#include "DataStructures/TempBuffer.hpp"
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/Christoffel.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/SpacetimeDerivativeOfSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalOneForm.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeNormalVector.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpatialMetric.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

// IWYU pragma: no_forward_declare Tensor

namespace GeneralizedHarmonic {
template <size_t SpatialDim, typename Frame, typename DataType>
void spacetime_riemann_tensor(
    const gsl::not_null<tnsr::abcc<DataType, VolumeDim, Frame>*> riemann,
    const tnsr::aa<DataType, SpatialDim, Frame> pi,
    const tnsr::iaa<DataType, VolumeDim, Frame>& phi,
    const tnsr::iaa<DataType, VolumeDim, Frame>& deriv_pi,
    const tnsr::ijaa<DataType, VolumeDim, Frame>& deriv_phi,
    const tnsr::II<DataType, VolumeDim, Frame>& inverse_spatial_metric) {
//   gr::spatial_metric(spatial_metric, spacetime_metric);
//   determinant_and_inverse(det_spatial_metric, inverse_spatial_metric,
//                           *spatial_metric);
//   gr::shift(shift, spacetime_metric, *inverse_spatial_metric);
//   gr::lapse(lapse, *shift, spacetime_metric);
//   gr::inverse_spacetime_metric(inverse_spacetime_metric, *lapse, *shift,
//                                *inverse_spatial_metric);
//   GeneralizedHarmonic::spacetime_derivative_of_spacetime_metric(
//       da_spacetime_metric, *lapse, *shift, pi, phi);
//   gr::christoffel_first_kind(christoffel_first_kind, *da_spacetime_metric);
//   raise_or_lower_first_index(christoffel_second_kind,
//                              *christoffel_first_kind,
//                              *inverse_spacetime_metric);
//   trace_last_indices(trace_christoffel, *christoffel_first_kind,
//                      *inverse_spacetime_metric);
//   gr::spacetime_normal_vector(normal_spacetime_vector, *lapse, *shift);
//   gr::spacetime_normal_one_form(normal_spacetime_one_form, *lapse);
}

template <size_t SpatialDim, typename Frame, typename DataType>
tnsr::abcc<DataType, VolumeDim, Frame> spacetime_riemann_tensor(
    const tnsr::aa<DataType, SpatialDim, Frame> pi,
    const tnsr::iaa<DataType, VolumeDim, Frame>& phi,
    const tnsr::iaa<DataType, VolumeDim, Frame>& deriv_pi,
    const tnsr::ijaa<DataType, VolumeDim, Frame>& deriv_phi,
    const tnsr::II<DataType, VolumeDim, Frame>& inverse_spatial_metric) {
  //
}
}  // namespace GeneralizedHarmonic
