// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"

/// \cond
namespace domain {
class BlockId;
}  // namespace domain
class DataVector;
template <size_t VolumeDim>
class Domain;
template <typename IdType, typename DataType>
class IdPair;
/// \endcond

/// \ingroup ComputationalDomainGroup
///
/// Computes the block logical coordinates and the containing `BlockId` of
/// a set of points, given coordinates in a particular frame.
///
/// \details Returns a std::vector<std::optional<IdPair<BlockId,coords>>>,
/// where the vector runs over the points and is indexed in the same order as
/// the input coordinates `x`. For each point, the `IdPair` holds the
/// block logical coords of that point and the `BlockId` of the `Block` that
/// contains that point.
/// The std::optional is invalid if the point is not in any Block.
/// If a point is on a shared boundary of two or more `Block`s, it is
/// returned only once, and is considered to belong to the `Block`
/// with the smaller `BlockId`.
///
/// \warning Since map inverses can involve numerical roundoff error, care must
/// be taken with points on shared block boundaries. They will be assigned to
/// the first block (by block ID) that contains the point _within roundoff
/// error_. Therefore, be advised to use the logical coordinates returned by
/// this function, which are guaranteed to be in [-1, 1] and can be safely
/// passed along to `element_logical_coordinates`.
template <size_t Dim, typename Frame>
auto block_logical_coordinates(
    const Domain<Dim>& domain, const tnsr::I<DataVector, Dim, Frame>& x,
    double time = std::numeric_limits<double>::signaling_NaN(),
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time = std::unordered_map<
            std::string,
            std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>{})
    -> std::vector<std::optional<
        IdPair<domain::BlockId, tnsr::I<double, Dim, ::Frame::BlockLogical>>>>;
