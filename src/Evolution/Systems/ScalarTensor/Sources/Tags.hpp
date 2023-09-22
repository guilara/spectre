// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/String.hpp"

namespace ScalarTensor::Sources {
namespace Tags {
struct ScalarSource : db::SimpleTag {
  using type = Scalar<DataVector>;
};

// Extra compute tags for debugging
/*!
 * \brief The \f$ f'(\Psi) \f$ term.
 */
struct CouplingFunctionDerivative : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "CouplingFunctionDerivative"; }
};

/*!
 * \brief The \f$ 8 (E^2 - B^2) \f$ term.
 */
struct GBScalar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "GBScalar"; }
};

/*!
 * \brief Rhs of the Psi equation.
 */
struct RhsPsi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "Rhs(Psi(CurvedScalarWave))"; }
};

/*!
 * \brief Rhs of the Pi equation.
 */
struct RhsPi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "Rhs(Pi(CurvedScalarWave))"; }
};

/*!
 * \brief Rhs of the Pi equation.
 */
struct RhsPhi : db::SimpleTag {
  using type = tnsr::i<DataVector, 3>;
  static std::string name() { return "Rhs(Phi(CurvedScalarWave))"; }
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

/*!
 * \brief Coupling parameter 1.
 */
struct ScalarFirstCouplingParameter {
  static std::string name() { return "EtaCouplingParameter"; }
  using type = double;
  static constexpr Options::String help{
      "First coupling parameter entering in the scalar field source, in code "
      "units"};
};

/*!
 * \brief Coupling parameter 2.
 */
struct ScalarSecondCouplingParameter {
  static std::string name() { return "ZetaCouplingParameter"; }
  using type = double;
  static constexpr Options::String help{
      "Second coupling parameter entering in the scalar field source, in code "
      "units"};
};
}  // namespace OptionTags

namespace Tags {
struct ScalarMass : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarMass>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double mass_psi) { return mass_psi; }
};

struct ScalarFirstCouplingParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarFirstCouplingParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double first_coupling_psi) {
    return first_coupling_psi;
  }
};

struct ScalarSecondCouplingParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarSecondCouplingParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double second_coupling_psi) {
    return second_coupling_psi;
  }
};
}  // namespace Tags
}  // namespace ScalarTensor::Sources
