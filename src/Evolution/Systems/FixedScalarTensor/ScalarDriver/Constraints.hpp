// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Constraints.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::ScalarDriver {

namespace Tags {
/*!
 * \brief Compute item to get the one-index constraint for the scalar driver
 * evolution system.
 *
 * \details See `one_index_constraint()`. Can be retrieved using
 * `CurvedScalarWave::Tags::OneIndexConstraint`.
 */
struct OneIndexConstraintCompute : OneIndexConstraint, db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<Psi, tmpl::size_t<3_st>, Frame::Inertial>,
                 Phi<3_st>>;
  using return_type = tnsr::i<DataVector, 3_st, Frame::Inertial>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>&) =
      &CurvedScalarWave::one_index_constraint<3_st>;
  using base = OneIndexConstraint;
};

/*!
 * \brief Compute item to get the two-index constraint for the scalar driver
 * evolution system.
 *
 * \details See `two_index_constraint()`. Can be retrieved using
 * `CurvedScalarWave::Tags::TwoIndexConstraint`.
 */
struct TwoIndexConstraintCompute : TwoIndexConstraint, db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<Phi<3_st>, tmpl::size_t<3_st>, Frame::Inertial>>;
  using return_type = tnsr::ij<DataVector, 3_st, Frame::Inertial>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result,
      const tnsr::ij<DataVector, 3_st, Frame::Inertial>&) =
      &CurvedScalarWave::two_index_constraint<3_st>;
  using base = TwoIndexConstraint;
};
}  // namespace Tags

}  // namespace fe::ScalarDriver
