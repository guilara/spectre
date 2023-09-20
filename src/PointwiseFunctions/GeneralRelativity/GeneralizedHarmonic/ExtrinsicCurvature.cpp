// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/ExtrinsicCurvature.hpp"

#include "DataStructures/DataVector.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

// IWYU pragma: no_forward_declare Tensor

namespace gh {
template <typename DataType, size_t SpatialDim, typename Frame>
void extrinsic_curvature(
    const gsl::not_null<tnsr::ii<DataType, SpatialDim, Frame>*> ex_curv,
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::aa<DataType, SpatialDim, Frame>& pi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi) {
  destructive_resize_components(ex_curv,
                                get_size(get<0>(spacetime_normal_vector)));
  for (auto& component : *ex_curv) {
    component = 0.0;
  }
  for (size_t i = 0; i < SpatialDim; ++i) {
    for (size_t j = i; j < SpatialDim; ++j) {
      for (size_t a = 0; a <= SpatialDim; ++a) {
        ex_curv->get(i, j) += 0.5 *
                              (phi.get(i, j + 1, a) + phi.get(j, i + 1, a)) *
                              spacetime_normal_vector.get(a);
      }
      ex_curv->get(i, j) += 0.5 * pi.get(i + 1, j + 1);
    }
  }
}

template <typename DataType, size_t SpatialDim, typename Frame>
tnsr::ii<DataType, SpatialDim, Frame> extrinsic_curvature(
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::aa<DataType, SpatialDim, Frame>& pi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi) {
  tnsr::ii<DataType, SpatialDim, Frame> ex_curv{};
  extrinsic_curvature(make_not_null(&ex_curv), spacetime_normal_vector, pi,
                      phi);
  return ex_curv;
}

template <typename DataType, size_t SpatialDim, typename Frame>
void grad_extrinsic_curvature(
    const gsl::not_null<tnsr::ijj<DataType, SpatialDim, Frame>*> grad_ex_curv,
    const tnsr::ijj<DataType, SpatialDim, Frame>& d_ex_curv,
    const tnsr::ii<DataType, SpatialDim, Frame>& ex_curv,
    const tnsr::Ijj<DataType, SpatialDim, Frame>&
        spatial_christoffel_second_kind) {
  destructive_resize_components(grad_ex_curv, get_size(get<0,0>(ex_curv)));
  tenex::evaluate<ti::i, ti::j, ti::k>(
      grad_ex_curv, d_ex_curv(ti::i, ti::j, ti::k) -
                        spatial_christoffel_second_kind(ti::L, ti::i, ti::j) *
                            ex_curv(ti::l, ti::k) -
                        spatial_christoffel_second_kind(ti::L, ti::i, ti::k) *
                            ex_curv(ti::j, ti::l));
}
}  // namespace gh

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(2, data)

#define INSTANTIATE(_, data)                                               \
  template void gh::extrinsic_curvature(                                   \
      const gsl::not_null<tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>*>  \
          ex_curv,                                                         \
      const tnsr::A<DTYPE(data), DIM(data), FRAME(data)>&                  \
          spacetime_normal_vector,                                         \
      const tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>& pi,             \
      const tnsr::iaa<DTYPE(data), DIM(data), FRAME(data)>& phi);          \
  template tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>                   \
  gh::extrinsic_curvature(                                                 \
      const tnsr::A<DTYPE(data), DIM(data), FRAME(data)>&                  \
          spacetime_normal_vector,                                         \
      const tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>& pi,             \
      const tnsr::iaa<DTYPE(data), DIM(data), FRAME(data)>& phi);          \
  template void gh::grad_extrinsic_curvature(                              \
      const gsl::not_null<tnsr::ijj<DTYPE(data), DIM(data), FRAME(data)>*> \
          grad_ex_curv,                                                    \
      const tnsr::ijj<DTYPE(data), DIM(data), FRAME(data)>& d_ex_curv,     \
      const tnsr::ii<DTYPE(data), DIM(data), FRAME(data)>& ex_curv,        \
      const tnsr::Ijj<DTYPE(data), DIM(data), FRAME(data)>&                \
          spatial_christoffel_second_kind);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (double, DataVector),
                        (Frame::Grid, Frame::Inertial))

#undef DIM
#undef DTYPE
#undef FRAME
#undef INSTANTIATE
