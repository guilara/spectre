// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"

namespace CurvedScalarWave::Sources {

namespace Tags {

struct ScalarSource : db::SimpleTag {
  using type = Scalar<DataVector>;
}

}  // namespace Tags

}  // namespace CurvedScalarWave::Sources
