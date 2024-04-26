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

// template <typename Frame>
void order_reduced_gb_scalar_with_tenex(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const tnsr::aa<DataVector, 3, Frame::Inertial>&
        trace_reversed_stress_energy,
    const tnsr::AA<DataVector, 3, Frame::Inertial>& inverse_spacetime_metric) {
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
                            inverse_spacetime_metric(ti::D, ti::C));
}

}  // namespace ScalarTensor
