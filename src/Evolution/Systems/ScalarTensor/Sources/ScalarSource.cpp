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
  scalar_source->get() += square(mass_psi) * get(psi);
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

void multiply_by_coupling_function_double_prime_quartic(
    gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi) {
  const auto ones_scalar = make_with_value<Scalar<DataVector>>(psi, 1.0);
  *scalar_source->get() *=
      -(first_coupling_psi / 4.0) * ones_scalar.get() -
      (3.0 * second_coupling_psi / 4.0) * square(psi.get());
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

void compute_rhs_psi(
    const gsl::not_null<Scalar<DataVector>*> dt_psi,
    const tnsr::i<DataVector, 3>& d_psi, const tnsr::i<DataVector, 3>& d_pi,
    const tnsr::ij<DataVector, 3>& d_phi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3>& phi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::i<DataVector, 3>& deriv_lapse,
    const tnsr::iJ<DataVector, 3>& deriv_shift,
    const tnsr::II<DataVector, 3>& upper_spatial_metric,
    const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2) {
  tenex::evaluate(dt_psi,
                  -lapse() * pi() + shift(ti::I) * d_psi(ti::i) +
                      gamma1() * shift(ti::J) * (d_psi(ti::j) - phi(ti::j)));
}

void compute_rhs_pi(
    const gsl::not_null<Scalar<DataVector>*> dt_pi,
    const tnsr::i<DataVector, 3>& d_psi, const tnsr::i<DataVector, 3>& d_pi,
    const tnsr::ij<DataVector, 3>& d_phi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3>& phi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::i<DataVector, 3>& deriv_lapse,
    const tnsr::iJ<DataVector, 3>& deriv_shift,
    const tnsr::II<DataVector, 3>& upper_spatial_metric,
    const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
    const Scalar<DataVector>& scalar_source) {
  tenex::evaluate(
      dt_pi,
      lapse() * pi() * trace_extrinsic_curvature() +
          shift(ti::I) * d_pi(ti::i) +
          lapse() * trace_spatial_christoffel(ti::I) * phi(ti::i) +
          gamma1() * gamma2() * shift(ti::I) * (d_psi(ti::i) - phi(ti::i)) -
          lapse() * upper_spatial_metric(ti::I, ti::J) * d_phi(ti::i, ti::j) -
          upper_spatial_metric(ti::I, ti::J) * phi(ti::i) * deriv_lapse(ti::j));
  // Add scalar source
  add_scalar_source_to_dt_pi_scalar(dt_pi, scalar_source, lapse);
}

void compute_rhs_phi(
    const gsl::not_null<tnsr::i<DataVector, 3, Frame::Inertial>*> dt_phi,
    const tnsr::i<DataVector, 3>& d_psi, const tnsr::i<DataVector, 3>& d_pi,
    const tnsr::ij<DataVector, 3>& d_phi, const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3>& phi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3>& shift,
    const tnsr::i<DataVector, 3>& deriv_lapse,
    const tnsr::iJ<DataVector, 3>& deriv_shift,
    const tnsr::II<DataVector, 3>& upper_spatial_metric,
    const tnsr::I<DataVector, 3>& trace_spatial_christoffel,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2) {
  tenex::evaluate<ti::i>(
      dt_phi, -lapse() * d_pi(ti::i) + shift(ti::J) * d_phi(ti::j, ti::i) +
                  gamma2() * lapse() * (d_psi(ti::i) - phi(ti::i)) -
                  pi() * deriv_lapse(ti::i) +
                  phi(ti::j) * deriv_shift(ti::i, ti::J));
}

void compute_order_reduced_scalar_curvature_source(
    gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& order_reduced_gb_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi) {
  // Make sure it has the same size
  // *scalar_source = make_with_value<Scalar<DataVector>>(psi, 0.);
  // Compute the Riemann squared scalar in vacuum
  scalar_source->get() = get(order_reduced_gb_scalar);
  // Multiply by the source coupling function
  multiply_by_coupling_function_prime_quartic(
      scalar_source, psi, first_coupling_psi, second_coupling_psi);

  // Add mass term
  scalar_source->get() += square(mass_psi) * psi.get();
}

}  // namespace ScalarTensor
