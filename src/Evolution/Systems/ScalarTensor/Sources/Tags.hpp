// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/String.hpp"

namespace ScalarTensor {
namespace Tags {
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
 * \brief The GB scalar term with nonvaccuum contributions from the scalar.
 */
struct OrderReducedGBScalar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "OrderReducedGBScalar"; }
};

}  // namespace Tags

namespace OptionTags {

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
}  // namespace ScalarTensor
