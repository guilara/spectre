// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines DataBox tags for scalar tensor system

#pragma once

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/TMPL.hpp"
namespace ScalarTensor {
namespace Tags {
/// Represents the trace reversed stress-energy tensor of the scalar
/// sector of the ScalarTensor system
struct TraceReversedStressEnergy : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
};

/// We use this prefix tag to avoid ambiguities with the ::gh system when
/// parsing observed variables
template <typename Tag>
struct CSW : db::PrefixTag, Tag {
  using type = typename Tag::type;
  using tag = Tag;  // remove_all_prefixes is called by Variables,
                    // SliceVariables, Subcell projection
                    // Is this a problem?
  static std::string name() { return "CSW(" + Tag::name() + ")"; }
};

} // namespace Tags
} // namespace ScalarTensor
