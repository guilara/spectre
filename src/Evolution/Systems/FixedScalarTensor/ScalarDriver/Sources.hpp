// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Options/Options.hpp"

namespace fe::ScalarDriver::Sources {

void add_scalar_driver_source_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi,
    const Scalar<DataVector>& scalar_driver_source,
    const Scalar<DataVector>& lapse);

void add_scalar_driver_friction_term_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi,
    const Scalar<DataVector>& scalar_driver_pi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3_st>& shift, const double scalar_tau_parameter,
    const double scalar_sigma_parameter);

void compute_scalar_driver_source(const gsl::not_null<return_type*> result,
                                  const Scalar<DataVector>& psi,
                                  const Scalar<DataVector>& target_psi);

void compute_target_psi(const gsl::not_null<return_type*> result);

namespace Tags {

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <size_t SpatialDim, typename Frame, typename DataType>
struct ScalarDriverSourceCompute : ScalarDriverSource, db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarDriver::Psi, fe::ScalarDriver::TargetPsi>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataVector>&,
      const Scalar<DataVector>&) = &compute_scalar_driver_source;
  using base = ScalarSource;
};

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <size_t SpatialDim, typename Frame, typename DataType>
struct TargetPsiCompute : TargetPsi, db::ComputeTag {
  using argument_tags = tmpl::list<fe::ScalarDriver::Psi>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>&) =
      &compute_target_psi;
  using base = ScalarSource;
};

}  // namespace Tags

}  // namespace fe::ScalarDriver::Sources
