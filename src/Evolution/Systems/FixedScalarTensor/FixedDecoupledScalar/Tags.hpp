// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Tags.hpp"

namespace fe::DecoupledScalar {
namespace Tags {}  // namespace Tags
namespace OptionTags {
/*!
 * \ingroup OptionGroupsGroup
 * Groups option tags related to the FixedDecoupledScalar evolution system.
 */
struct Group {
  static std::string name() { return "FixedDecoupledScalar"; }
  static constexpr Options::String help{
      "Options for the FixedDecoupledScalar evolution system"};
  using group = evolution::OptionTags::SystemGroup;
};
}  // namespace OptionTags
}  // namespace fe::DecoupledScalar
