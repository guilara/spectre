// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <optional>
#include <string>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/BoundaryCorrection.hpp"
#include "Evolution/Systems/ScalarTensor/BoundaryCorrections/ProductOfCorrections.hpp"
#include "Evolution/Systems/ScalarTensor/System.hpp"
#include "Evolution/Systems/ScalarTensor/TimeDerivativeTerms.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/BoundaryCorrections.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Options/Options.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace {
template <typename DerivedCorrection, typename PackagedFieldTagList,
          typename EvolvedTagList, typename FluxTagList, typename TempTagList,
          // Should we keep this?
          typename PrimTagList,
          //
          typename VolumeTagList>
struct ComputeBoundaryCorrectionHelperImpl;

template <typename DerivedCorrection, typename... PackagedFieldTags,
          typename... EvolvedTags, typename... FluxTags, typename... TempTags,
          // Should we keep this?
          typename... PrimTags,
          //
          typename... VolumeTags>
struct ComputeBoundaryCorrectionHelperImpl<
    DerivedCorrection, tmpl::list<PackagedFieldTags...>,
    tmpl::list<EvolvedTags...>, tmpl::list<FluxTags...>,
    tmpl::list<TempTags...>,
    // Should we keep this?
    tmpl::list<PrimTags...>,
    //
    tmpl::list<VolumeTags...>> {
  template <typename PackagedVariables, typename EvolvedVariables,
            typename FluxVariables, typename TempVariables,
            typename PrimVariables, typename VolumeVariables>
  static double
  dg_package_data(/* Add variables */
                  const gsl::not_null<PackagedVariables*> packaged_variables,
                  const EvolvedVariables& evolved_variables,
                  const FluxVariables& flux_variables,
                  const TempVariables& temp_variables,
                  // Should we keep this?
                  const PrimVariables& prim_variables,
                  //
                  const VolumeVariables& volume_variables,
                  const tnsr::i<DataVector, 3, Frame::Inertial>&
                      normal_covector,
                  const tnsr::I<DataVector, 3, Frame::Inertial>& normal_vector,
                  const std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>&
                      mesh_velocity,
                  const std::optional<Scalar<DataVector>>&
                      normal_dot_mesh_velocity,
                  const DerivedCorrection& derived_correction) {
    return derived_correction.dg_package_data(
        make_not_null(&get<PackagedFieldTags>(*packaged_variables))...,
        get<EvolvedTags>(evolved_variables)...,
        get<FluxTags>(flux_variables)..., get<TempTags>(temp_variables)...,
        // Should we keep this?
        get<PrimTags>(prim_variables)...,
        //
        normal_covector, normal_vector, mesh_velocity, normal_dot_mesh_velocity,
        get<VolumeTags>(volume_variables)...);
  }

  template <typename EvolvedVariables, typename PackagedVariables>
  static void dg_boundary_terms(
      const gsl::not_null<EvolvedVariables*> boundary_corrections,
      const PackagedVariables& internal_packaged_fields,
      const PackagedVariables& external_packaged_fields,
      dg::Formulation dg_formulation,
      const DerivedCorrection& derived_correction) {
    derived_correction.dg_boundary_terms(
        make_not_null(&get<EvolvedTags>(*boundary_corrections))...,
        get<PackagedFieldTags>(internal_packaged_fields)...,
        get<PackagedFieldTags>(external_packaged_fields)..., dg_formulation);
  }
};

template <typename DerivedCorrection, typename EvolvedTagList,
          typename FluxTagList>
using ComputeBoundaryCorrectionHelper = ComputeBoundaryCorrectionHelperImpl<
    DerivedCorrection, typename DerivedCorrection::dg_package_field_tags,
    EvolvedTagList, FluxTagList,
    typename DerivedCorrection::dg_package_data_temporary_tags,
    // Should we keep this?
    typename DerivedCorrection::dg_package_data_primitive_tags,
    //
    typename DerivedCorrection::dg_package_data_volume_tags>;

template <typename DerivedGhCorrection, typename DerivedScalarCorrection>
void test_boundary_correction_combination(
    const DerivedGhCorrection& derived_gh_correction,
    const DerivedScalarCorrection& derived_scalar_correction,
    const ScalarTensor::BoundaryCorrections::ProductOfCorrections<
        DerivedGhCorrection, DerivedScalarCorrection>&
        derived_product_correction,
    const dg::Formulation formulation) {
  CHECK(derived_product_correction.gh_correction() == derived_gh_correction);
  CHECK(derived_product_correction.scalar_correction() ==
        derived_scalar_correction);
  // ... much more
}

} // namespace

SPECTRE_TEST_CASE(
    "Unit.ScalarTensor.BoundaryCorrections.ProductOfCorrections",
    "[Unit][Evolution]") {
  // scoped to separate out each product combination
  {
    INFO("Product correction UpwindPenalty and UpwindPenalty");
    // Check that UpwindPenalty for CurvedScalarWave
    // does not need template parameters
    CurvedScalarWave::BoundaryCorrections::UpwindPenalty<3_st>
        scalar_correction{};
    GeneralizedHarmonic::BoundaryCorrections::UpwindPenalty<3_st>
        gh_correction{};
    TestHelpers::test_creation<std::unique_ptr<
        CurvedScalarWave::BoundaryCorrections::BoundaryCorrection>>(
        "ProductUpwindPenaltyAndUpwindPenalty:\n"
        "  UpwindPenalty:\n"
        "  UpwindPenalty:");
    ScalarTensor::BoundaryCorrections::ProductOfCorrections<
        GeneralizedHarmonic::BoundaryCorrections::UpwindPenalty<3_st>,
        CurvedScalarWave::BoundaryCorrections::UpwindPenalty<3_st>>
        product_boundary_correction{gh_correction, scalar_correction};
    for (const auto formulation :
         {dg::Formulation::StrongInertial, dg::Formulation::WeakInertial}) {
      test_boundary_correction_combination<
          GeneralizedHarmonic::BoundaryCorrections::UpwindPenalty<3_st>,
          CurvedScalarWave::BoundaryCorrections::UpwindPenalty<3_st>>(
          gh_correction, scalar_correction, product_boundary_correction,
          formulation);
    }
}
