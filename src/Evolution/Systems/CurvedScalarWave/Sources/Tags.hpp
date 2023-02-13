// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/Options.hpp"

namespace CurvedScalarWave::Sources {

namespace Tags {

struct ScalarSource : db::SimpleTag {
  using type = Scalar<DataVector>;
};

} // namespace Tags

namespace OptionTags {

/*!
 * \brief Scalar mass parameter.
 */
struct ScalarMass {
  static std::string name() { return "ScalarMass"; }
  using type = double;
  static constexpr Options::String help{
      "Mass of the scalar field in code units"};
  // using group = XXX;
};

}  // namespace OptionTags

namespace Tags {

struct ScalarMass : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarMass>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double mass_psi) { return mass_psi; }
};

} // namespace Tags

}  // namespace CurvedScalarWave::Sources
