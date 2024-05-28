// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Tags.hpp"

namespace fe::sgb {
namespace Tags {}  // namespace Tags
namespace OptionTags {
/*!
 * \ingroup OptionGroupsGroup
 * Groups option tags related to the Fixed SGB evolution system.
 */
struct Group {
  static std::string name() { return "FixedSGB"; }
  static constexpr Options::String help{
      "Options for the fixed SGB evolution system"};
  using group = evolution::OptionTags::SystemGroup;
};
}  // namespace OptionTags
}  // namespace fe::sgb
