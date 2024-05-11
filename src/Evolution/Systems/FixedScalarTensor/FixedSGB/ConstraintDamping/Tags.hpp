// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/ConstraintDamping/DampingFunction.hpp"
#include "Options/String.hpp"

/// \cond
namespace fe::sgb::OptionTags {
struct Group;
}  // namespace fe::sgb::OptionTags
/// \endcond

namespace fe::sgb::ConstraintDamping {
namespace OptionTags {

// template <size_t VolumeDim, typename Fr>
// struct DampingFunctionGamma1 {
//   using type = std::unique_ptr<
//       ::fe::sgb::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
//   static constexpr Options::String help{
//       "DampingFunction for damping parameter gamma1"};
//   using group = fe::sgb::OptionTags::Group;
// };

// template <size_t VolumeDim, typename Fr>
// struct DampingFunctionGamma2 {
//   using type = std::unique_ptr<
//       ::fe::sgb::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
//   static constexpr Options::String help{
//       "DampingFunction for damping parameter gamma2"};
//   using group = fe::sgb::OptionTags::Group;
// };

template <size_t VolumeDim, typename Fr>
struct DampingFunctionSigmaParameter {
  using type = std::unique_ptr<
      ::fe::sgb::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  static constexpr Options::String help{
      "DampingFunction for driver parameter sigma"};
  using group = fe::sgb::OptionTags::Group;
};

template <size_t VolumeDim, typename Fr>
struct DampingFunctionTauParameter {
  using type = std::unique_ptr<
      ::fe::sgb::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  static constexpr Options::String help{
      "DampingFunction for driver parameter tau"};
  using group = fe::sgb::OptionTags::Group;
};

}  // namespace OptionTags

namespace Tags {

/*!
 * \brief A DampingFunction to compute the scalar driver parameter
 * \f$\sigma\f$.
 */
template <size_t VolumeDim, typename Fr>
struct DampingFunctionSigmaParameter : db::SimpleTag {
  using DampingFunctionType =
      ::fe::sgb::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::fe::sgb::ConstraintDamping::OptionTags::
                     DampingFunctionSigmaParameter<VolumeDim, Fr>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& damping_function) {
    return damping_function->get_clone();
  }
};

/*!
 * \brief A DampingFunction to compute the scalar driver parameter
 * \f$\tau\f$.
 */
template <size_t VolumeDim, typename Fr>
struct DampingFunctionTauParameter : db::SimpleTag {
  using DampingFunctionType =
      ::fe::sgb::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::fe::sgb::ConstraintDamping::OptionTags::
                     DampingFunctionTauParameter<VolumeDim, Fr>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& damping_function) {
    return damping_function->get_clone();
  }
};

}  // namespace Tags
}  // namespace fe::sgb::ConstraintDamping
