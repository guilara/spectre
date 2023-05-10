// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/Options.hpp"

namespace ScalarTensor::Sources {

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

/*!
 * \brief Constraint damping parameter corresponding to gamma2.
 */
struct ConstraintDampingExternalParameterGamma2 {
  static std::string name() {
    return "ConstraintDampingExternalParameterGamma2";
  }
  using type = double;
  static constexpr Options::String help{
      "Constraint damping parameter corresponding to gamma2"};
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

struct ConstraintDampingExternalParameterGamma2 : db::SimpleTag {
  using type = double;
  using option_tags =
      tmpl::list<OptionTags::ConstraintDampingExternalParameterGamma2>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double external_gamma2_parameter) {
    return external_gamma2_parameter;
  }
};

} // namespace Tags

}  // namespace ScalarTensor::Sources
