// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Sources/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Sources {

/*!
 * \brief Compute the scalar source term for the CurvedScalarWave system.
 *
 * \details The scalar source term depends on the problem at hand. For a scalar
 * subject to a potential,
 *
 * \f[
 * \mathrm{scalar source} = \partial V / \partial \psi .
 * \f]
 *
 * For a massive scalar \f$ \partial V / \partial \psi = m_{\psi}^2 \psi \f$.
 */
void compute_scalar_source(gsl::not_null<Scalar<DataVector>*> scalar_source,
                           const Scalar<DataVector>& psi);

/*!
 * \brief Add in the scalar source term for the CurvedScalarWave system.
 *
 * \details The only stress energy source term in the CurvedScalarWave
 * evolution equations is in the equation for \f$\Pi\f$:
 *
 * \f[
 * \partial_t \Pi_{ab} = \alpha (\mathrm{scalar source})
 * + \text{(other CurvedScalarWave source terms)} ,
 * \f]
 *
 * where \f$ \mathrm{scalar source} = \partial V / \partial \psi \f$.
 *
 * For a massive scalar \f$ \partial V / \partial \psi = m_{\psi}^2 \psi \f$.
 *
 * This function adds that contribution to the existing value of `dt_pi`. The
 * spacetime terms in the Curved Scalar equation should be computed before
 * passing the `dt_pi` to this function for updating.
 *
 * \see `CurvedScalarWave::TimeDerivative` for details about the spacetime
 * part of the time derivative calculation.
 */
void add_scalar_source_to_dt_pi(gsl::not_null<Scalar<DataVector>*> dt_pi,
                                const Scalar<DataVector>& scalar_source,
                                const Scalar<DataVector>& lapse);

namespace Tags {

/*!
 * \brief Compute tag for the scalar source.
 *
 * \details Call compute_scalar_source.
 */
struct ScalarSourceCompute : ScalarSource, db::ComputeTag {
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Psi>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>&) =
      &compute_scalar_source;
  using base = ScalarSource;
};

}  // namespace Tags

}  // namespace CurvedScalarWave::Sources
