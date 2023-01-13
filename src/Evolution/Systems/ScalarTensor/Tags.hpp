// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
// Add scalar tags?
#include "Utilities/TMPL.hpp"

namespace ScalarTensor {

namespace Tags {
/// %Tags for the combined system of the Generalized Harmonic formulation for
/// the Einstein field equations and the scalar field system.

namespace detail {
// A const reference to another tag, used for rerouting arguments in the
// combined system utilities
template <typename Tag, typename Type = db::const_item_type<Tag, tmpl::list<>>>
struct TemporaryReference {
  using tag = Tag;
  using type = const Type&;
};
}  // namespace detail

/// Represents the trace reversed stress-energy tensor of the scalar
/// sector of the ScalarTensor system
struct TraceReversedStressEnergy : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
};

/// Represents the stress-energy tensor of the scalar of the
/// ScalarTensor system
struct StressEnergy : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
};

} // namespace Tags
} // namespace ScalarTensor
