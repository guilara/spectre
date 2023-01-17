// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <type_traits>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/TempBuffer.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TypeTraits/IsA.hpp"

/// @{
/*!
 * \ingroup DataStructuresGroup
 * \brief Retrieves a desired tag from data structures containing tags.
 *
 * \details Given multiple containers retrieve a desired tag from the first
 * container in which it is found. The containers are searched in the order
 * in which they are supplied (i.e. `first_vars` is checked before `vars`).
 * An error will be emitted if the tag cannot be found in any container.
 */
template <typename Tag, typename... TagsOrTagsList, typename... Ts,
          template <class...> class U,
          Requires<(not tt::is_a_v<gsl::not_null, U<TagsOrTagsList...>>)and not(
              ... and tt::is_a_v<gsl::not_null, std::decay_t<Ts>>)> = nullptr>
const typename Tag::type& get(const U<TagsOrTagsList...>& first_vars,
                              const Ts&... vars) {
  if constexpr (tt::is_a_v<db::DataBox, U<TagsOrTagsList...>>) {
    if constexpr (not db::tag_is_retrievable_v<Tag, U<TagsOrTagsList...>> and
                  sizeof...(Ts) != 0) {
      return get<Tag>(vars...);
    } else {
      return get<Tag>(first_vars);
    }
  } else {
    if constexpr (not tmpl::list_contains_v<
                      tmpl::flatten<tmpl::list<TagsOrTagsList...>>, Tag> and
                  sizeof...(Ts) != 0) {
      return get<Tag>(vars...);
    } else {
      return get<Tag>(first_vars);
    }
  }
}

template <typename Tag, typename... TagsOrTagsList, typename... Ts,
          template <class...> class U>
gsl::not_null<typename Tag::type*> get(
    const gsl::not_null<U<TagsOrTagsList...>*> first_vars,
    const gsl::not_null<Ts>... vars) {
  if constexpr (tt::is_a_v<db::DataBox, U<TagsOrTagsList...>>) {
    if constexpr (not db::tag_is_retrievable_v<Tag, U<TagsOrTagsList...>> and
                  sizeof...(Ts) != 0) {
      if constexpr (sizeof...(Ts) == 1) {
        return make_not_null(&get<Tag>(*vars...));
      } else {
        return get<Tag>(vars...);
      }
    } else {
      return make_not_null(&get<Tag>(*first_vars));
    }
  } else {
    if constexpr (not tmpl::list_contains_v<
                      tmpl::flatten<tmpl::list<TagsOrTagsList...>>, Tag> and
                  sizeof...(Ts) != 0) {
      (void)first_vars;
      if constexpr (sizeof...(Ts) == 1) {
        return make_not_null(&get<Tag>(*vars...));
      } else {
        return get<Tag>(vars...);
      }
    } else {
      return make_not_null(&get<Tag>(*first_vars));
    }
  }
}
/// @}
