// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/SetNumberOfGridPoints.hpp"

namespace ScalarTensor {

void add_scalar_source_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
    const Scalar<DataVector>& scalar_source, const Scalar<DataVector>& lapse) {
  get(*dt_pi_scalar) += get(lapse) * get(scalar_source);
}

void mass_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& psi, const double mass_psi) {
  get(*scalar_source) = square(mass_psi) * get(psi);
}

void compute_scalar_curvature_source(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi) {
  // Make sure it has the same size
  // *scalar_source = make_with_value<Scalar<DataVector>>(psi, 0.);
  // Compute the Riemann squared scalar in vacuum
  scalar_source->get() =
      8.0 * (weyl_electric_scalar.get() - weyl_magnetic_scalar.get());
  // Multiply by the source coupling function
  multiply_by_coupling_function_prime_quartic(
      scalar_source, psi, first_coupling_psi, second_coupling_psi);

  // Add mass term
  scalar_source->get() += square(mass_psi) * psi.get();
}

Scalar<DataVector> compute_scalar_curvature_source(
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi) {
  Scalar<DataVector> result{};
  compute_scalar_curvature_source(make_not_null(&result), weyl_electric_scalar,
                                  weyl_magnetic_scalar, psi, first_coupling_psi,
                                  second_coupling_psi, mass_psi);
  return result;
}

void multiply_by_coupling_function_prime_quartic(
    gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi) {
  *scalar_source->get() *= -(first_coupling_psi / 4.0) * psi.get() -
                           (second_coupling_psi / 4.0) * cube(psi.get());
}

// Extra functions for debugging
void compute_coupling_function_derivative(
    const gsl::not_null<Scalar<DataVector>*> result,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi) {
  *result = make_with_value<Scalar<DataVector>>(psi, 1.0);
  multiply_by_coupling_function_prime_quartic(result, psi, first_coupling_psi,
                                              second_coupling_psi);
}

void compute_gb_scalar(const gsl::not_null<Scalar<DataVector>*> gb_scalar,
                       const Scalar<DataVector>& weyl_electric_scalar,
                       const Scalar<DataVector>& weyl_magnetic_scalar) {
  // *gb_scalar = make_with_value<Scalar<DataVector>>(weyl_electric_scalar, 0.);
  // Compute the Riemann squared scalar in vacuum
  gb_scalar->get() =
      8.0 * (weyl_electric_scalar.get() - weyl_magnetic_scalar.get());
}

}  // namespace ScalarTensor
