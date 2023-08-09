// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

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
namespace ScalarTensor {
namespace Tags {
/*!
 * \brief Represents the trace-reversed stress-energy tensor of the scalar
 * field.
 */
template <typename DataType, size_t Dim, typename Fr = Frame::Inertial>
struct TraceReversedStressEnergy : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Fr>;
};

/*!
 * \brief Tag holding the source term of the scalar equation.
 *
 * \details This tag hold the source term \f$ \mathcal{S} \f$,
 * entering a wave equation of the form
 * \f[
 *   \Box \Psi = \mathcal{S} ~.
 * \f]
 */
struct ScalarSource : db::SimpleTag {
  using type = Scalar<DataVector>;
};

}  // namespace Tags

namespace OptionTags {
/*!
 * \brief Scalar mass parameter.
 */
struct ScalarMass {
  static std::string name() { return "ScalarMass"; }
  using type = double;
  static constexpr Options::String help{
      "Mass of the scalar field in code units"};
};
}  // namespace OptionTags

namespace Tags {
struct ScalarMass : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarMass>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double mass_psi) { return mass_psi; }
};

/*!
 * \brief Prefix tag to avoid ambiguities when observing variables with the same
 * name in both parent systems.
 * \note Since we also add compute tags for these quantities, we do not make
 * this a derived class of `Tag`. Otherwise, we would have tags with repeated
 * base tags in the `ObservationBox`.
 */
template <typename Tag>
struct CSW : db::PrefixTag, db::SimpleTag {
  using type = typename Tag::type;
  using tag = Tag;
};

/*!
 * \brief Compute tag for observing the evolved variable Psi.
 * \details This tag copies the data from the evolved variable Tag
 * CurvedScalarWave::Tags::Psi into the wrapped tag used for observation.
 */
struct CSWPsiCompute : CSW<CurvedScalarWave::Tags::Psi>, db::ComputeTag {
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Psi>;
  using return_type = Scalar<DataVector>;
  static constexpr void function(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>& psi) {
    get(*result) = get(psi);
  }
  using base = CSW<CurvedScalarWave::Tags::Psi>;
};

/*!
 * \brief Compute tag for observing the evolved variable Pi.
 * \details This tag copies the data from the evolved variable Tag
 * CurvedScalarWave::Tags::Pi into the wrapped tag used for observation.
 */
struct CSWPiCompute : CSW<CurvedScalarWave::Tags::Pi>, db::ComputeTag {
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Pi>;
  using return_type = Scalar<DataVector>;
  static constexpr void function(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>& pi) {
    get(*result) = get(pi);
  }
  using base = CSW<CurvedScalarWave::Tags::Pi>;
};

/*!
 * \brief Compute tag for observing the evolved variable Phi.
 * \details This tag copies the data from the evolved variable Tag
 * CurvedScalarWave::Tags::Phi into the wrapped tag used for observation.
 */
template <size_t Dim>
struct CSWPhiCompute : CSW<CurvedScalarWave::Tags::Phi<Dim>>, db::ComputeTag {
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Phi<Dim>>;
  using return_type = tnsr::i<DataVector, Dim>;
  static constexpr void function(const gsl::not_null<return_type*> result,
                                    const tnsr::i<DataVector, Dim>& phi) {
    for (size_t i = 0; i < Dim; ++i) {
      result->get(i) = phi.get(i);
    }
  }
  using base = CSW<CurvedScalarWave::Tags::Phi<Dim>>;
};

/*!
 * \brief Computes the scalar-wave one-index constraint.
 * \details The one-index constraint is assigned to a wrapped tag to avoid
 * clashes with the ::gh constraints during observation.
 */
template <size_t Dim>
struct CSWOneIndexConstraintCompute
    : CSW<CurvedScalarWave::Tags::OneIndexConstraint<Dim>>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<CurvedScalarWave::Tags::Psi, tmpl::size_t<Dim>,
                               Frame::Inertial>,
                 CurvedScalarWave::Tags::Phi<Dim>>;
  using return_type = tnsr::i<DataVector, Dim>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const tnsr::i<DataVector, Dim>&,
                                    const tnsr::i<DataVector, Dim>&) =
      &CurvedScalarWave::one_index_constraint<Dim>;
  using base = CSW<CurvedScalarWave::Tags::OneIndexConstraint<Dim>>;
};

/*!
 * \brief Computes the scalar-wave two-index constraint.
 * \details The two-index constraint is assigned to a wrapped tag to avoid
 * clashes with the ::gh constraints during observation.
 */
template <size_t Dim>
struct CSWTwoIndexConstraintCompute
    : CSW<CurvedScalarWave::Tags::TwoIndexConstraint<Dim>>,
      db::ComputeTag {
  using argument_tags =
      tmpl::list<::Tags::deriv<CurvedScalarWave::Tags::Phi<Dim>,
                               tmpl::size_t<Dim>, Frame::Inertial>>;
  using return_type = tnsr::ij<DataVector, Dim>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const tnsr::ij<DataVector, Dim>&) =
      &CurvedScalarWave::two_index_constraint<Dim>;
  using base = CSW<CurvedScalarWave::Tags::TwoIndexConstraint<Dim>>;
};

}  // namespace Tags

}  // namespace ScalarTensor
