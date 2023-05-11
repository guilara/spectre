// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor {

/*!
 * \brief Add in the trace-reversed stress-energy source term to the \f$\Pi\f$
 * evolved variable of the Generalized Harmonic system.
 *
 * \details The only stress energy source term in the Generalized Harmonic
 * evolution equations is in the equation for \f$\Pi_{a b}\f$:
 *
 * \f[
 * \partial_t \Pi_{ab} + \text{(spatial derivative terms)} =
 * \text{(GH source terms)}
 * - 16 \pi \alpha (T_{ab} - \frac{1}{2} g_{a b} T^c{}_c)
 * \f]
 *
 * (note that this function takes as argument the trace-reversed stress energy
 * tensor)
 *
 * This function adds that contribution to the existing value of `dt_pi`. The
 * spacetime terms in the GH equation should be computed before passing the
 * `dt_pi` to this function for updating.
 *
 * \see `GeneralizedHarmonic::TimeDerivative` for details about the spacetime
 * part of the time derivative calculation.
 */
void add_stress_energy_term_to_dt_pi(
    gsl::not_null<tnsr::aa<DataVector, 3>*> dt_pi,
    const tnsr::aa<DataVector, 3>& trace_reversed_stress_energy,
    const Scalar<DataVector>& lapse);

/*!
 * \brief Calculate the trace-reversed stress-energy tensor \f$(T_{\mu \nu} -
 * 1/2 g_{\mu \nu} g^{\lambda \sigma} T_{\lambda \sigma}) \f$ associated with
 * the scalar part of the ScalarTensor system.
 *
 * \details The stress energy tensor is needed to compute the backreaction of
 * the scalar to the spacetime degrees of freedom. The stress energy is
 * set to zero for now.
 *
 * \f[
 * T_{\mu \nu} = 0,
 * \f]
 *
 * where ...
 */
void trace_reversed_stress_energy(
    gsl::not_null<tnsr::aa<DataVector, 3>*> stress_energy,
    /* Add scalar and scalar gradients */
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3>& phi_scalar,
    const tnsr::aa<DataVector, 3, Frame::Inertial>& spacetime_metric,
    // const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
    const Scalar<DataVector>& lapse);
// void trace_reversed_stress_energy(
//     gsl::not_null<tnsr::aa<DataVector, 3>*> stress_energy);

namespace Tags {

/*!
 * \brief Compute tag for the trace reversed stress energy tensor.
 *
 * \details Call trace_reversed_stress_energy.
 */
struct TraceReversedStressEnergyCompute : TraceReversedStressEnergy,
                                          db::ComputeTag {
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Pi, CurvedScalarWave::Tags::Phi<3_st>,
                 gr::Tags::SpacetimeMetric<3_st, ::Frame::Inertial, DataVector>,
                 gr::Tags::Lapse<DataVector>>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, 3>*> result,
      /* Add scalar and scalar gradients */
      const Scalar<DataVector>&, const tnsr::i<DataVector, 3>&,
      const tnsr::aa<DataVector, 3, Frame::Inertial>&,
      const Scalar<DataVector>&) = &trace_reversed_stress_energy;
  using base = TraceReversedStressEnergy;
};
}  // namespace Tags

} // namespace ScalarTensor
