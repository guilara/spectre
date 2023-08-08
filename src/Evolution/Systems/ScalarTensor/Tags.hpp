// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
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

/// We use this prefix tag to avoid ambiguities with the ::gh system when
/// parsing observed variables
template <typename Tag>
struct CSW : db::PrefixTag, Tag {
  using type = typename Tag::type;
  using tag = Tag;
  static std::string name() { return "CSW(" + Tag::name() + ")"; }
};

/*!
 * \brief Compute tag for observing the evolved variable Psi.
 */
struct CSWPsiCompute : CSW<CurvedScalarWave::Tags::Psi>, db::ComputeTag {
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Psi>;
  using return_type = Scalar<DataVector>;
  static constexpr void function(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>& input) {
    get(*result) = get(input);
  }
  using base = CSW<CurvedScalarWave::Tags::Psi>;
};

/*!
 * \brief Compute tag for observing the evolved variable Pi.
 */
struct CSWPiCompute : CSW<CurvedScalarWave::Tags::Pi>, db::ComputeTag {
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Pi>;
  using return_type = Scalar<DataVector>;
  static constexpr void function(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>& input) {
    get(*result) = get(input);
  }
  using base = CSW<CurvedScalarWave::Tags::Pi>;
};

/*!
 * \brief Compute tag for observing the evolved variable Phi.
 */
template <size_t Dim>
struct CSWPhiCompute : CSW<CurvedScalarWave::Tags::Phi<Dim>>, db::ComputeTag {
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Phi<Dim>>;
  using return_type = tnsr::i<DataVector, Dim>;
  static constexpr void function(const gsl::not_null<return_type*> result,
                                    const tnsr::i<DataVector, Dim>& input) {
    for (size_t i = 0; i < Dim; ++i) {
      result->get(i) = input.get(i);
    }
  }
  using base = CSW<CurvedScalarWave::Tags::Phi<Dim>>;
};

}  // namespace Tags

}  // namespace ScalarTensor
