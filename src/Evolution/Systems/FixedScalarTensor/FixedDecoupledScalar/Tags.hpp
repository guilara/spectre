// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once
#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Tags.hpp"
#include "Options/String.hpp"

/// \cond
class DataVector;
template <class>
class Variables;
/// \endcond

namespace fe::DecoupledScalar {
namespace Tags {

struct ReSource : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct ImSource : db::SimpleTag {
  using type = Scalar<DataVector>;
};

}  // namespace Tags
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
