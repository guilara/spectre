// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/FixedSGB/Constraints.hpp"

#include <cstddef>

#include "DataStructures/LeviCivitaIterator.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Constraints.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::sgb {
template <typename DataType, size_t SpatialDim, typename Frame>
void f_constraint(
    const gsl::not_null<tnsr::a<DataType, SpatialDim, Frame>*> constraint,
    const tnsr::a<DataType, SpatialDim, Frame>& gauge_function,
    const tnsr::ab<DataType, SpatialDim, Frame>& spacetime_d_gauge_function,
    const tnsr::a<DataType, SpatialDim, Frame>& spacetime_normal_one_form,
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::II<DataType, SpatialDim, Frame>& inverse_spatial_metric,
    const tnsr::AA<DataType, SpatialDim, Frame>& inverse_spacetime_metric,
    const tnsr::aa<DataType, SpatialDim, Frame>& pi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& phi,
    const tnsr::iaa<DataType, SpatialDim, Frame>& d_pi,
    const tnsr::ijaa<DataType, SpatialDim, Frame>& d_phi,
    const Scalar<DataType>& gamma2,
    const tnsr::iaa<DataType, SpatialDim, Frame>& three_index_constraint,
    const tnsr::aa<DataType, SpatialDim, Frame>& trace_reversed_stress_energy,
    const tnsr::aa<DataType, SpatialDim, Frame>& tensor_driver) {
  gh::f_constraint<DataType, SpatialDim, Frame>(
      constraint, gauge_function, spacetime_d_gauge_function,
      spacetime_normal_one_form, spacetime_normal_vector,
      inverse_spatial_metric, inverse_spacetime_metric, pi, phi, d_pi, d_phi,
      gamma2, three_index_constraint);
  fe::sgb::f_constraint_add_stress_energy_term(
      constraint, inverse_spacetime_metric, spacetime_normal_vector,
      spacetime_normal_one_form, trace_reversed_stress_energy, tensor_driver);
}

template <typename DataType, size_t SpatialDim, typename Frame>
void f_constraint_add_stress_energy_term(
    const gsl::not_null<tnsr::a<DataType, SpatialDim, Frame>*> constraint,
    const tnsr::AA<DataType, SpatialDim, Frame>& inverse_spacetime_metric,
    const tnsr::A<DataType, SpatialDim, Frame>& spacetime_normal_vector,
    const tnsr::a<DataType, SpatialDim, Frame>& spacetime_normal_one_form,
    const tnsr::aa<DataType, SpatialDim, Frame>& trace_reversed_stress_energy,
    const tnsr::aa<DataType, SpatialDim, Frame>& tensor_driver) {
  // This term, like many terms in the f constraint, may benefit from
  // allocating a temporary for the trace. However, once we apply that
  // optimization it should be applied to all terms that can benefit from
  // temporary storage, and the allocation shared among all of the terms.

  // Define kappa = 8.0 * M_PI * G
  const double kappa = 1.0;

  for (size_t a = 0; a < SpatialDim + 1; ++a) {
    for (size_t b = 0; b < SpatialDim + 1; ++b) {
      constraint->get(a) -= 2.0 * kappa * spacetime_normal_vector.get(b) *
                            trace_reversed_stress_energy.get(a, b);
      // Tensor driver term
      constraint->get(a) -= 2.0 * kappa * spacetime_normal_vector.get(b) *
                            tensor_driver.get(a, b);
      for (size_t c = 0; c < SpatialDim + 1; ++c) {
        constraint->get(a) += kappa * spacetime_normal_one_form.get(a) *
                              inverse_spacetime_metric.get(b, c) *
                              trace_reversed_stress_energy.get(b, c);
        // Tensor driver term
        constraint->get(a) += kappa * spacetime_normal_one_form.get(a) *
                              inverse_spacetime_metric.get(b, c) *
                              tensor_driver.get(b, c);
      }
    }
  }
}

}  // namespace fe::sgb

// Explicit Instantiations
#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(2, data)

#define INSTANTIATE(_, data)                                              \
  template void fe::sgb::f_constraint(                                    \
      const gsl::not_null<tnsr::a<DTYPE(data), DIM(data), FRAME(data)>*>  \
          constraint,                                                     \
      const tnsr::a<DTYPE(data), DIM(data), FRAME(data)>& gauge_function, \
      const tnsr::ab<DTYPE(data), DIM(data), FRAME(data)>&                \
          spacetime_d_gauge_function,                                     \
      const tnsr::a<DTYPE(data), DIM(data), FRAME(data)>&                 \
          spacetime_normal_one_form,                                      \
      const tnsr::A<DTYPE(data), DIM(data), FRAME(data)>&                 \
          spacetime_normal_vector,                                        \
      const tnsr::II<DTYPE(data), DIM(data), FRAME(data)>&                \
          inverse_spatial_metric,                                         \
      const tnsr::AA<DTYPE(data), DIM(data), FRAME(data)>&                \
          inverse_spacetime_metric,                                       \
      const tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>& pi,            \
      const tnsr::iaa<DTYPE(data), DIM(data), FRAME(data)>& phi,          \
      const tnsr::iaa<DTYPE(data), DIM(data), FRAME(data)>& d_pi,         \
      const tnsr::ijaa<DTYPE(data), DIM(data), FRAME(data)>& d_phi,       \
      const Scalar<DTYPE(data)>& gamma2,                                  \
      const tnsr::iaa<DTYPE(data), DIM(data), FRAME(data)>&               \
          three_index_constraint,                                         \
      const tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>&                \
          trace_reversed_stress_energy,                                   \
      const tnsr::aa<DTYPE(data), DIM(data), FRAME(data)>& tensor_driver);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (double, DataVector),
                        (Frame::Grid, Frame::Inertial))

#undef DIM
#undef DTYPE
#undef FRAME
#undef INSTANTIATE
