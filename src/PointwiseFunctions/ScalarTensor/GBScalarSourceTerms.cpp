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
void gb_scalar(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector> weyl_electric_scalar,
    const Scalar<DataVector> weyl_magnetic_scalar,
    const tnsr::aa<DataVector, 3, Frame> trace_reversed_stress_energy,
    const tnsr::AA<DataType, 3, Frame> inverse_spacetime_metric) {
  static constexpr double kappa = 8.0 * M_PI;
  static constexpr double two_over_three = 2.0 / 3.0;
  get(*result) = 8.0 * (weyl_electric_scalar - weyl_magnetic_scalar);
  // Raise index of the trace reversed stress energy tensor
  tnsr::Aa<DataVector, 3, Frame> trace_reversed_stress_energy_up_down =
      make_with_value<tnsr::Aa<DataVector, 3, Frame>>(
          get<0, 0>(inverse_spacetime_metric), 0.0);
  for (size_t b = 0; b < 4; ++b) {
    for (size_t c = 0; c < 4; ++c) {
      for (size_t a = 0; a < 4; ++a) {
        trace_reversed_stress_energy_up_down.get(b, c) +=
            inverse_spacetime_metric.get(b, a) *
            trace_reversed_stress_energy.get(a, c);
      }
    }
  }
  // Compute the trace of the trace reversed stress energy tensor
  Scalar<DataVector> trace_of_trace_reversed_stress_energy =
      make_with_value<Scalar<DataVector>>(get<0, 0>(inverse_spacetime_metric),
                                          0.0);
  for (size_t a = 0; a < 4; ++a) {
    get(trace_of_trace_reversed_stress_energy) +=
        trace_reversed_stress_energy_up_down.get(a, a);
  }
  // Compute the full contraction of the trace reversed stress energy tensor
  // with itself
  for (size_t a = 0; a < 4; ++a) {
    for (size_t b = 0; b < 4; ++b) {
      get(*result) += -2.0 * square(kappa) *
                      trace_reversed_stress_energy_up_down.get(a, b) *
                      trace_reversed_stress_energy_up_down.get(b, a);
    }
  }
  // Add square of trace of trace reversed stress energy tensor
  get(*result) += two_over_three * square(kappa) *
                  square(trace_of_trace_reversed_stress_energy);
}

template <typename Frame>
void gb_scalar_with_tenex(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector> weyl_electric_scalar,
    const Scalar<DataVector> weyl_magnetic_scalar,
    const tnsr::aa<DataVector, 3, Frame> trace_reversed_stress_energy,
    const tnsr::AA<DataType, 3, Frame> inverse_spacetime_metric) {
  static constexpr double kappa = 8.0 * M_PI;
  static constexpr double two_over_three = 2.0 / 3.0;
  tenex::evaluate(result,
                  // Weyl squared in terms of electric and magnetic scalars
                  8.0 * weyl_electric_scalar() -
                      8.0 * weyl_magnetic_scalar()
                      // Trace reversed stress energy squared
                      - 2.0 * square(kappa) *
                            trace_reversed_stress_energy(ti::a, ti::b) *
                            inverse_spacetime_metric(ti::B, ti::C) *
                            trace_reversed_stress_energy(ti::c, ti::d) *
                            inverse_spacetime_metric(ti::D, ti::A)
                      // Square of the trace of the trace reversed stress energy
                      + two_over_three * square(kappa) *
                            trace_reversed_stress_energy(ti::a, ti::b) *
                            inverse_spacetime_metric(ti::B, ti::A) *
                            trace_reversed_stress_energy(ti::c, ti::d) *
                            inverse_spacetime_metric(ti::D, ti::C))

}  // namespace ScalarTensor
