// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/ConstraintDamping/DampingFunction.hpp"
#include "Options/String.hpp"

/// \cond
namespace fe::DecoupledScalar::OptionTags {
struct Group;
}  // namespace fe::DecoupledScalar::OptionTags
/// \endcond

namespace fe::DecoupledScalar::ConstraintDamping {
namespace OptionTags {

template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma1 {
  using type = std::unique_ptr<
      ::fe::DecoupledScalar::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  static constexpr Options::String help{
      "DampingFunction for damping parameter gamma1"};
  using group = fe::DecoupledScalar::OptionTags::Group;
};

template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma2 {
  using type = std::unique_ptr<
      ::fe::DecoupledScalar::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  static constexpr Options::String help{
      "DampingFunction for damping parameter gamma2"};
  using group = fe::DecoupledScalar::OptionTags::Group;
};

template <size_t VolumeDim, typename Fr>
struct DampingFunctionScalarSigmaParameter {
  using type = std::unique_ptr<
      ::fe::DecoupledScalar::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  static constexpr Options::String help{
      "DampingFunction for driver parameter sigma"};
  using group = fe::DecoupledScalar::OptionTags::Group;
};

template <size_t VolumeDim, typename Fr>
struct DampingFunctionScalarTauParameter {
  using type = std::unique_ptr<
      ::fe::DecoupledScalar::ConstraintDamping::DampingFunction<VolumeDim, Fr>>;
  static constexpr Options::String help{
      "DampingFunction for driver parameter tau"};
  using group = fe::DecoupledScalar::OptionTags::Group;
};

}  // namespace OptionTags

namespace Tags {

/*!
 * \brief A DampingFunction to compute the constraint damping parameter
 * \f$\gamma_0\f$.
 */
template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma1 : db::SimpleTag {
  using DampingFunctionType =
      ::fe::DecoupledScalar::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::fe::DecoupledScalar::ConstraintDamping::OptionTags::
                     DampingFunctionGamma1<VolumeDim, Fr>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& damping_function) {
    return damping_function->get_clone();
  }
};

/*!
 * \brief A DampingFunction to compute the constraint damping parameter
 * \f$\gamma_0\f$.
 */
template <size_t VolumeDim, typename Fr>
struct DampingFunctionGamma2 : db::SimpleTag {
  using DampingFunctionType =
      ::fe::DecoupledScalar::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::fe::DecoupledScalar::ConstraintDamping::OptionTags::
                     DampingFunctionGamma2<VolumeDim, Fr>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& damping_function) {
    return damping_function->get_clone();
  }
};

/*!
 * \brief A DampingFunction to compute the scalar driver parameter
 * \f$\sigma\f$.
 */
template <size_t VolumeDim, typename Fr>
struct DampingFunctionScalarSigmaParameter : db::SimpleTag {
  using DampingFunctionType =
      ::fe::DecoupledScalar::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::fe::DecoupledScalar::ConstraintDamping::OptionTags::
                     DampingFunctionScalarSigmaParameter<VolumeDim, Fr>>;

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
struct DampingFunctionScalarTauParameter : db::SimpleTag {
  using DampingFunctionType =
      ::fe::DecoupledScalar::ConstraintDamping::DampingFunction<VolumeDim, Fr>;
  using type = std::unique_ptr<DampingFunctionType>;
  using option_tags =
      tmpl::list<::fe::DecoupledScalar::ConstraintDamping::OptionTags::
                     DampingFunctionScalarTauParameter<VolumeDim, Fr>>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& damping_function) {
    return damping_function->get_clone();
  }
};

}  // namespace Tags
}  // namespace fe::DecoupledScalar::ConstraintDamping
