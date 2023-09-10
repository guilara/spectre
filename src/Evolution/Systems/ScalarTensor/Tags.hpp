// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines DataBox tags for scalar tensor system

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/TMPL.hpp"
#include "Options/String.hpp"


namespace ScalarTensor {
namespace OptionTags {

/*!
 * \brief Amplitude of the gaussian function for the constraint damping
 * parameter for the scalar driver equation.
 */
struct AmplitudeConstraintGamma2 {
  static std::string name() {
    return "AmplitudeConstraintGamma2(CurvedScalarWave)";
  }
  using type = double;
  static constexpr Options::String help{
      "Amplitude of the gaussian function for the constraint damping parameter "
      "for the scalar driver equation."};
  static type lower_bound() { return 0.; }
};

/*!
 * \brief Width of the gaussian function for the constraint damping
 * parameter for the scalar driver equation.
 */
struct SigmaConstraintGamma2 {
  static std::string name() {
    return "SigmaConstraintGamma2(CurvedScalarWave)";
  }
  using type = double;
  static constexpr Options::String help{
      "Width of the gaussian function for the constraint damping parameter "
      "for the scalar driver equation."};
  static type lower_bound() { return 0.; }
};

/*!
 * \brief Asymptotic value of the gaussian function for the constraint damping
 * parameter for the scalar driver equation.
 */
struct OffsetConstraintGamma2 {
  static std::string name() {
    return "OffsetConstraintGamma2(CurvedScalarWave)";
  }
  using type = double;
  static constexpr Options::String help{
      "Asymptotic value for the gaussian function for the constraint damping "
      "parameter for the scalar driver equation."};
  static type lower_bound() { return 0.; }
};

}  // namespace OptionTags

namespace Tags {
/// Represents the trace reversed stress-energy tensor of the scalar
/// sector of the ScalarTensor system
struct TraceReversedStressEnergy : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
};

struct AmplitudeConstraintGamma2 : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::AmplitudeConstraintGamma2>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double amplitude_gaussian) {
    return amplitude_gaussian;
  }
};

struct SigmaConstraintGamma2 : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::SigmaConstraintGamma2>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double sigma_gaussian) {
    return sigma_gaussian;
  }
};

struct OffsetConstraintGamma2 : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::OffsetConstraintGamma2>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double offset_gaussian) {
    return offset_gaussian;
  }
};

} // namespace Tags
} // namespace ScalarTensor
