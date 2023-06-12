// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/Options.hpp"

namespace fe::ScalarDriver::Sources {
namespace Tags {

struct ScalarDriverSource : db::SimpleTag {
  using type = Scalar<DataVector>;
};

}  // namespace Tags

namespace OptionTags {
/*!
 * \brief Scalar sigma parameter.
 */
struct ScalarSigmaParameter {
  static std::string name() { return "ScalarSigmaParameter"; }
  using type = double;
  static constexpr Options::String help{
      "Sigma parameter for the scalar diver in code units"};
};

/*!
 * \brief Scalar tau parameter.
 */
struct ScalarTauParameter {
  static std::string name() { return "ScalarTauParameter"; }
  using type = double;
  static constexpr Options::String help{
      "Tau parameter for the scalar driver in code units"};
};
}  // namespace OptionTags

namespace Tags {
struct ScalarSigmaParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarSigmaParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double scalar_sigma_parameter) {
    return scalar_sigma_parameter;
  }
};

struct ScalarTauParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarTauParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double scalar_tau_parameter) {
    return scalar_tau_parameter;
  }
};

}  // namespace Tags
}  // namespace fe::ScalarDriver::Sources
