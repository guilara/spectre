// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Options/Options.hpp"

namespace fe::ScalarDriver::Sources {

void compute_scalar_driver_source(const gsl::not_null<return_type*> result,
                                  const Scalar<DataVector>& psi,
                                  const Scalar<DataVector>& target_psi);

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

}  // namespace Tags

}  // namespace fe::ScalarDriver::Sources
