// Distributed under the MIT License.
// See LICENSE.txt for details.

///\file
/// Defines functions to calculate the scalar wave constraints in
/// curved spacetime

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Constraints.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarTensor {
namespace Tags {
/*!
 * \brief Compute item to get the one-index constraint for the scalar-wave
 * evolution system.
 *
 * \details See `one_index_constraint()`. Can be retrieved using
 * `CurvedScalarWave::Tags::OneIndexConstraint`.
 */
template <size_t SpatialDim>
struct OneIndexConstraintCompute
    : CSW<CurvedScalarWave::Tags::OneIndexConstraint<SpatialDim>>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<CSW<CurvedScalarWave::Tags::Psi>,
                               tmpl::size_t<SpatialDim>, Frame::Inertial>,
                 CSW<CurvedScalarWave::Tags::Phi<SpatialDim>>>;
  using return_type = tnsr::i<DataVector, SpatialDim, Frame::Inertial>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result,
      const tnsr::i<DataVector, SpatialDim, Frame::Inertial>&,
      const tnsr::i<DataVector, SpatialDim, Frame::Inertial>&) =
      &CurvedScalarWave::one_index_constraint<SpatialDim>;
  using base = CSW<OneIndexConstraint<SpatialDim>>;
};

/*!
 * \brief Compute item to get the two-index constraint for the scalar-wave
 * evolution system.
 *
 * \details See `two_index_constraint()`. Can be retrieved using
 * `CurvedScalarWave::Tags::TwoIndexConstraint`.
 */
template <size_t SpatialDim>
struct TwoIndexConstraintCompute
    : CSW<CurvedScalarWave::Tags::TwoIndexConstraint<SpatialDim>>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<CSW<CurvedScalarWave::Tags::Phi<SpatialDim>>,
                               tmpl::size_t<SpatialDim>, Frame::Inertial>>;
  using return_type = tnsr::ij<DataVector, SpatialDim, Frame::Inertial>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result,
      const tnsr::ij<DataVector, SpatialDim, Frame::Inertial>&) =
      &CurvedScalarWave::two_index_constraint<SpatialDim>;
  using base = CSW<TwoIndexConstraint<SpatialDim>>;
};
}  // namespace Tags
}  // namespace ScalarTensor
