// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <pup.h>

#include "Evolution/Systems/GeneralizedHarmonic/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/System.hpp"
#include "Evolution/Systems/CurvedScalarWave/BoundaryCorrections/Factory.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarTensor::BoundaryCorrections {

namespace detail {

template <typename DerivedGhCorrection, typename DerivedScalarCorrection,
          typename GhPackageFieldTagList, typename ScalarPackageFieldTagList,
          typename GhEvolvedTagList, typename ScalarEvolvedTags,
          typename GhFluxTagList, typename ScalarFluxTagList,
          typename GhTempTagList, typename ScalarTempTagList,
          typename DeduplicatedTempTags, typename GhPrimTagList,
          typename ScalarPrimTagList, typename GhVolumeTagList,
          typename ScalarVolumeTagList>
struct ProductOfCorrectionsImpl;

template <typename DerivedGhCorrection, typename DerivedScalarCorrection,
          typename... GhPackageFieldTags, typename... ScalarPackageFieldTags,
          typename... GhEvolvedTags, typename... ScalarEvolvedTags,
          typename... GhFluxTags, typename... ScalarFluxTags,
          typename... GhTempTags, typename... ScalarTempTags,
          // What are the following tags ?
            typename... DeduplicatedTempTags, typename... GhPrimTags,
            typename... ScalarPrimTags,
          //
          typename... GhVolumeTags, typename... ScalarVolumeTags>
struct ProductOfCorrectionsImpl<
    DerivedGhCorrection, DerivedScalarCorrection,
    tmpl::list<GhPackageFieldTags...>, tmpl::list<ScalarPackageFieldTags...>,
    tmpl::list<GhEvolvedTags...>, tmpl::list<ScalarEvolvedTags...>,
    tmpl::list<GhFluxTags...>, tmpl::list<ScalarFluxTags...>,
    tmpl::list<GhTempTags...>, tmpl::list<ScalarTempTags...>,
    // What are the following tags ?
    tmpl::list<DeduplicatedTempTags...>, tmpl::list<GhPrimTags...>,
    tmpl::list<ScalarPrimTags...>,
    //
    tmpl::list<GhVolumeTags...>, tmpl::list<ScalarVolumeTags...>> {
  static double dg_package_data(
                  const gsl::not_null<
                      typename GhPackageFieldTags::type*>... gh_packaged_fields,
                  const gsl::not_null<typename ScalarPackageFieldTags::
                                          type*>... scalar_packaged_fields,
                  const typename GhEvolvedTags::type&... gh_variables,
                  const typename ScalarEvolvedTags::type&... scalar_variables,
                  const typename GhFluxTags::type&... gh_fluxes,
                  const typename ScalarFluxTags::type&... scalar_fluxes,
                  // What are the following tags ?
                  const typename DeduplicatedTempTags::type&... temporaries,
                  const typename GhPrimTags::type&... gh_primitives,
                  const typename ScalarPrimTags::type&... scalar_primitives,
                  //
                  const tnsr::i<DataVector, 3, Frame::Inertial>&
                      normal_covector,
                  const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,
                  const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
                      mesh_velocity,
                  const std::optional<Scalar<DataVector>>&
                      normal_dot_mesh_velocity,
                  const typename GhVolumeTags::type&... gh_volume_quantities,
                  const typename ScalarVolumeTags::type&...
                      scalar_volume_quantities,
                  const DerivedGhCorrection& gh_correction,
                  const DerivedScalarCorrection& scalar_correction){
      tuples::TaggedTuple<
          Tags::detail::TemporaryReference<DeduplicatedTempTags>...>
          shuffle_refs{temporaries...};
      return std::max(
          gh_correction.dg_package_data(
              gh_packaged_fields..., gh_variables..., gh_fluxes...,
              tuples::get<Tags::detail::TemporaryReference<GhTempTags>>(
                  shuffle_refs)...,
              gh_primitives..., normal_covector, normal_vector, mesh_velocity,
              normal_dot_mesh_velocity, gh_volume_quantities...),
          scalar_correction.dg_package_data(
              scalar_packaged_fields..., scalar_variables...,
              scalar_fluxes...,
              tuples::get<Tags::detail::TemporaryReference<ScalarTempTags>>(
                  shuffle_refs)...,
              scalar_primitives..., normal_covector, normal_vector,
              mesh_velocity, normal_dot_mesh_velocity,
              scalar_volume_quantities...));
  }

  static void dg_boundary_terms(
      const gsl::not_null<
          typename GhEvolvedTags::type*>... gh_boundary_corrections,
      const gsl::not_null<
          typename ScalarEvolvedTags::type*>... scalar_boundary_corrections,
      const typename GhPackageFieldTags::type&... gh_internal_packaged_fields,
      const typename ScalarPackageFieldTags::type&...
            scalar_internal_packaged_fields,
      const typename GhPackageFieldTags::type&... gh_external_packaged_fields,
      const typename ScalarPackageFieldTags::type&...
            scalar_external_packaged_fields,
      const dg::Formulation dg_formulation,
      const DerivedGhCorrection& gh_correction,
      const DerivedScalarCorrection& scalar_correction){
        /* Define */
      gh_correction.dg_boundary_terms(
          gh_boundary_corrections..., gh_internal_packaged_fields...,
          gh_external_packaged_fields..., dg_formulation);
      scalar_correction.dg_boundary_terms(
          scalar_boundary_corrections...,
          scalar_internal_packaged_fields...,
          scalar_external_packaged_fields..., dg_formulation);
  }
};
}  // namespace detail

/*!
 * \brief Apply a boundary condition to the combined Generalized Harmonic (GH)
 * and scalar field system using boundary corrections defined separately.
 */
template <typename DerivedGhCorrection, typename DerivedScalarCorrection>
class ProductOfCorrections final : public BoundaryCorrection {
 public:
  // Define
  using dg_package_field_tags =
      tmpl::append<typename DerivedGhCorrection::dg_package_field_tags,
                   typename DerivedScalarCorrection::dg_package_field_tags>;

  using dg_package_data_temporary_tags = tmpl::remove_duplicates<tmpl::append<
      typename DerivedGhCorrection::dg_package_data_temporary_tags,
      typename DerivedScalarCorrection::dg_package_data_temporary_tags>>;

//   using dg_package_data_primitive_tags =
//       typename DerivedScalarCorrection::dg_package_data_primitive_tags;

  using dg_package_data_volume_tags = tmpl::append<
      typename DerivedGhCorrection::dg_package_data_volume_tags,
      typename DerivedScalarCorrection::dg_package_data_volume_tags>;

  using derived_product_correction_impl = detail::ProductOfCorrectionsImpl<
      DerivedGhCorrection, DerivedScalarCorrection,
      typename DerivedGhCorrection::dg_package_field_tags,
      typename DerivedScalarCorrection::dg_package_field_tags,
      typename GeneralizedHarmonic::System<3_st>::variables_tag::tags_list,
      typename CurvedScalarWave::System<3_st>::variables_tag::tags_list,
      db::wrap_tags_in<
          ::Tags::Flux,
          typename GeneralizedHarmonic::System<3_st>::flux_variables,
          tmpl::size_t<3_st>, Frame::Inertial>,
      db::wrap_tags_in<::Tags::Flux,
                       typename CurvedScalarWave::System<3_st>::flux_variables,
                       tmpl::size_t<3_st>, Frame::Inertial>,
      typename DerivedGhCorrection::dg_package_data_temporary_tags,
      typename DerivedScalarCorrection::dg_package_data_temporary_tags,
      dg_package_data_temporary_tags, tmpl::list<>,
    //   typename DerivedScalarCorrection::dg_package_data_primitive_tags,
      typename DerivedGhCorrection::dg_package_data_volume_tags,
      typename DerivedScalarCorrection::dg_package_data_volume_tags>;

  static std::string name() {
    return "Product" + pretty_type::name<DerivedGhCorrection>() + "And" +
           pretty_type::name<DerivedScalarCorrection>();
  }

  struct GhCorrection {
    using type = DerivedGhCorrection;
    static std::string name() {
      return pretty_type::name<DerivedGhCorrection>();
    }
    static constexpr Options::String help{
        "The Generalized Harmonic part of the product boundary condition"};
  };
  struct ScalarCorrection {
    using type = DerivedScalarCorrection;
    static std::string name() {
      return pretty_type::name<DerivedScalarCorrection>();
    }
    static constexpr Options::String help{
          "The Curved Scalar part of the product boundary condition"};
    };

  using options = tmpl::list<GhCorrection, ScalarCorrection>;

  static constexpr Options::String help = {
      "Direct product of a GH and Curved Scalar boundary correction. "
      "See the documentation for the two individual boundary corrections for "
      "further details."};

  ProductOfCorrections() = default;
  ProductOfCorrections(DerivedGhCorrection gh_correction,
                       DerivedScalarCorrection scalar_correction)
      : derived_gh_correction_{gh_correction},
        derived_scalar_correction_{scalar_correction} {}
  ProductOfCorrections(const ProductOfCorrections&) = default;
  ProductOfCorrections& operator=(const ProductOfCorrections&) = default;
  ProductOfCorrections(ProductOfCorrections&&) = default;
  ProductOfCorrections& operator=(ProductOfCorrections&&) = default;
  ~ProductOfCorrections() override = default;

  /// \cond
  explicit ProductOfCorrections(CkMigrateMessage* msg)
      : BoundaryCorrection(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ProductOfCorrections);  // NOLINT
  /// \endcond
  void pup(PUP::er& p) override {
    BoundaryCorrection::pup(p);
    p | derived_gh_correction_;
    p | derived_scalar_correction_;
  }

  std::unique_ptr<BoundaryCorrection> get_clone() const override {
    return std::make_unique<ProductOfCorrections>(*this);
  }

  template <typename... Args>
  double dg_package_data(Args&&... args) const {
    return derived_product_correction_impl::dg_package_data(
        std::forward<Args>(args)..., derived_gh_correction_,
        derived_scalar_correction_);
  }

  template <typename... Args>
  void dg_boundary_terms(Args&&... args) const {
    derived_product_correction_impl::dg_boundary_terms(
        std::forward<Args>(args)..., derived_gh_correction_,
        derived_scalar_correction_);
  }

  const DerivedGhCorrection& gh_correction() const {
    return derived_gh_correction_;
  }

  const DerivedScalarCorrection& scalar_correction() const {
    return derived_scalar_correction_;
  }

 private:
  DerivedGhCorrection derived_gh_correction_;
  DerivedScalarCorrection derived_scalar_correction_;
};

/// \cond
template <typename DerivedGhCorrection, typename DerivedScalarCorrection>
PUP::able::PUP_ID ProductOfCorrections<DerivedGhCorrection,
                                       DerivedScalarCorrection>::my_PUP_ID =
    0;  // NOLINT
/// \endcond
}  // namespace ScalarTensor::BoundaryCorrections
