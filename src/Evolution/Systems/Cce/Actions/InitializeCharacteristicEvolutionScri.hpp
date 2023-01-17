// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "Evolution/Systems/Cce/ScriPlusInterpolationManager.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/Rational.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

namespace Cce {

/// \cond
template <class Metavariables>
struct AnalyticWorldtubeBoundary;
/// \endcond
namespace Actions {

/*!
 * \ingroup ActionsGroup
 * \brief Initializes the `CharacteristicEvolution` component with contents
 * needed to perform the interpolation at scri+.
 *
 * \details Sets up the \ref DataBoxGroup to be ready to store data in the scri+
 * interpolators and perform interpolation for the final scri+ outputs.
 *
 * \ref DataBoxGroup changes:
 * - Modifies: nothing
 * - Adds:
 *  - `Cce::Tags::InterpolationManager<ComplexDataVector, Tag>` for each `Tag`
 * in `scri_values_to_observe`
 * - Removes: nothing
 */
template <typename ScriValuesToObserve, typename BoundaryComponent>
struct InitializeCharacteristicEvolutionScri {
  using simple_tags_from_options = tmpl::flatten<tmpl::list<
      InitializationTags::ScriInterpolationOrder,
      tmpl::conditional_t<
          tt::is_a_v<AnalyticWorldtubeBoundary, BoundaryComponent>,
          tmpl::list<Tags::AnalyticBoundaryDataManager>, tmpl::list<>>>>;

  using const_global_cache_tags =
      tmpl::list<Tags::LMax, Tags::NumberOfRadialPoints>;

  using simple_tags =
      tmpl::transform<ScriValuesToObserve,
                      tmpl::bind<Tags::InterpolationManager,
                                 tmpl::pin<ComplexDataVector>, tmpl::_1>>;

  using compute_tags = tmpl::list<>;

  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTags>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    initialize_impl(make_not_null(&box),
                    typename Metavariables::scri_values_to_observe{});
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }

  template <typename TagList, typename... TagPack>
  static void initialize_impl(const gsl::not_null<db::DataBox<TagList>*> box,
                              tmpl::list<TagPack...> /*meta*/) {
    const size_t target_number_of_points =
        db::get<InitializationTags::ScriInterpolationOrder>(*box);
    const size_t vector_size =
        Spectral::Swsh::number_of_swsh_collocation_points(
            db::get<Spectral::Swsh::Tags::LMaxBase>(*box));
    // silence compiler warnings when pack is empty
    (void)vector_size;
    if constexpr (sizeof...(TagPack) > 0) {
      Initialization::mutate_assign<simple_tags>(
          box, ScriPlusInterpolationManager<ComplexDataVector, TagPack>{
                   target_number_of_points, vector_size,
                   std::make_unique<intrp::BarycentricRationalSpanInterpolator>(
                       2 * target_number_of_points - 1,
                       2 * target_number_of_points + 2)}...);
    }
  }
};
}  // namespace Actions
}  // namespace Cce
