// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Constraints.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/String.hpp"
#include "Utilities/TMPL.hpp"

/*!
 * \brief Tags for the scalar tensor system.
 */
namespace ScalarTensor::TwoScalars {
namespace Tags {}  // namespace Tags

namespace OptionTags {}  // namespace OptionTags

namespace Tags {

/*!
 * \brief Prefix tag to avoid ambiguities when observing variables with the same
 * name all component systems.
 * \note Since we also add compute tags for these quantities, we do not make
 * this a derived class of `Tag`. Otherwise, we would have tags with repeated
 * base tags in the `ObservationBox`.
 */
template <typename Tag, size_t ScalarLabel>
struct Csw : db::PrefixTag, db::SimpleTag {
  using type = typename Tag::type;
  using tag = Tag;
};

/*!
 * \brief Computes the scalar-wave one-index constraint.
 * \details The one-index constraint is assigned to a wrapped tag to avoid
 * clashes with the ::gh constraints during observation.
 * \note We do not use ScalarTensor::Tags::CswCompute to retrieve the
 * CurvedScalarWave constraints since we do not want to add the bare compute
 * tags (::CurvedScalarWave::Tags::OneIndexConstraintCompute and
 * ::CurvedScalarWave::Tags::TwoIndexConstraintCompute) directly in the
 * ObservationBox, since that will make them observable and would lead to a
 * clash with the ::gh constraint tags.
 */
template <size_t ScalarLabel>
struct CswOneIndexConstraintCompute
    : Csw<CurvedScalarWave::Tags::OneIndexConstraint<3>>,
      db::ComputeTag {
  static constexpr size_t Dim = 3;
  using argument_tags =
      tmpl::list<::Tags::deriv<Csw<CurvedScalarWave::Tags::Psi, ScalarLabel>,
                               tmpl::size_t<Dim>, Frame::Inertial>,
                 Csw<CurvedScalarWave::Tags::Phi<Dim>, ScalarLabel>>;
  using return_type = tnsr::i<DataVector, Dim>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const tnsr::i<DataVector, Dim>&,
                                    const tnsr::i<DataVector, Dim>&) =
      &CurvedScalarWave::one_index_constraint<Dim>;
  using base =
      Csw<CurvedScalarWave::Tags::OneIndexConstraint<Dim>, ScalarLabel>;
};

/*!
 * \brief Computes the scalar-wave two-index constraint.
 * \details The two-index constraint is assigned to a wrapped tag to avoid
 * clashes with the ::gh constraints during observation.
 * \note We do not use ScalarTensor::Tags::CswCompute to retrieve the
 * CurvedScalarWave constraints since we do not want to add the bare compute
 * tags (::CurvedScalarWave::Tags::OneIndexConstraintCompute and
 * ::CurvedScalarWave::Tags::TwoIndexConstraintCompute) directly in the
 * ObservationBox, since that will make them observable and would lead to a
 * clash with the ::gh constraint tags.
 */
template <size_t ScalarLabel>
struct CswTwoIndexConstraintCompute
    : Csw<CurvedScalarWave::Tags::TwoIndexConstraint<3>, ScalarLabel>,
      db::ComputeTag {
  static constexpr size_t Dim = 3;
  using argument_tags = tmpl::list<
      ::Tags::deriv<Csw<CurvedScalarWave::Tags::Phi<Dim>, ScalarLabel>,
                    tmpl::size_t<Dim>, Frame::Inertial>>;
  using return_type = tnsr::ij<DataVector, Dim>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const tnsr::ij<DataVector, Dim>&) =
      &CurvedScalarWave::two_index_constraint<Dim>;
  using base =
      Csw<CurvedScalarWave::Tags::TwoIndexConstraint<Dim>, ScalarLabel>;
};

}  // namespace Tags

}  // namespace ScalarTensor::TwoScalars
