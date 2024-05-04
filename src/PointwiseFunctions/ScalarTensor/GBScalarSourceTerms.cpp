// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/ScalarTensor/GBScalarSourceTerms.hpp"

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Math.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarTensor {

template <typename Frame>
void order_reduced_gb_scalar_with_tenex(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const tnsr::aa<DataVector, 3, Frame>& trace_reversed_stress_energy,
    const tnsr::AA<DataVector, 3, Frame>& inverse_spacetime_metric) {
  static constexpr double two_over_three = 2.0 / 3.0;
  tenex::evaluate(result,
                  // Weyl squared in terms of electric and magnetic scalars
                  // with complement
                  8.0 * (weyl_electric_scalar() - weyl_magnetic_scalar())
                      // Trace reversed stress energy squared
                      - 2.0 * trace_reversed_stress_energy(ti::a, ti::b) *
                            inverse_spacetime_metric(ti::B, ti::C) *
                            trace_reversed_stress_energy(ti::c, ti::d) *
                            inverse_spacetime_metric(ti::D, ti::A)
                      // Square of the trace of the trace reversed stress energy
                      + two_over_three *
                            trace_reversed_stress_energy(ti::a, ti::b) *
                            inverse_spacetime_metric(ti::B, ti::A) *
                            trace_reversed_stress_energy(ti::c, ti::d) *
                            inverse_spacetime_metric(ti::D, ti::C));
}

}  // namespace ScalarTensor

#define FRAME(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                      \
  template void ScalarTensor::order_reduced_gb_scalar_with_tenex( \
      const gsl::not_null<Scalar<DataVector>*> result,            \
      const Scalar<DataVector>& weyl_electric_scalar,             \
      const Scalar<DataVector>& weyl_magnetic_scalar,             \
      const tnsr::aa<DataVector, 3, FRAME(data)>&                 \
          trace_reversed_stress_energy,                           \
      const tnsr::AA<DataVector, 3, FRAME(data)>& inverse_spacetime_metric);

GENERATE_INSTANTIATIONS(INSTANTIATE, (Frame::Grid, Frame::Inertial))

#undef FRAME
#undef INSTANTIATE
