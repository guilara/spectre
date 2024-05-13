// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/FaceNormal.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "ParallelAlgorithms/ApparentHorizonFinder/HorizonAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
namespace gsl {
template <class T>
class not_null;
}  // namespace gsl
namespace Tags {
template <typename Tag>
struct Normalized;
}  // namespace Tags
/// \endcond

// IWYU pragma: no_forward_declare Tensor

namespace fe::ScalarTensorDriver {
/// @{
/*!
 * \brief Compute the characteristic speeds for the scalar tensor driver system.
 */
template <size_t Dim, typename Frame>
void characteristic_speeds(
    const gsl::not_null<tnsr::a<DataVector, 3, Frame>*> char_speeds,
    // const Scalar<DataVector>& gamma_1,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame>& shift,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form,
    const std::optional<tnsr::I<DataVector, Dim, Frame>>& mesh_velocity);

template <size_t Dim, typename Frame>
std::array<DataVector, 4> characteristic_speeds(
    // const Scalar<DataVector>& gamma_1,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame>& shift,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form,
    const std::optional<tnsr::I<DataVector, Dim, Frame>>& mesh_velocity);

template <size_t Dim, typename Frame>
void characteristic_speeds(
    const gsl::not_null<std::array<DataVector, 4>*> char_speeds,
    // const Scalar<DataVector>& gamma_1,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame>& shift,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form,
    const std::optional<tnsr::I<DataVector, Dim, Frame>>& mesh_velocity);

template <size_t Dim, typename Frame>
struct CharacteristicSpeedsCompute
    : Tags::CharacteristicSpeeds<DataVector, Dim, Frame>,
      db::ComputeTag {
  using base = Tags::CharacteristicSpeeds<DataVector, Dim, Frame>;
  using type = typename base::type;
  using argument_tags = tmpl::list<
      //   ::gh::ConstraintDamping::Tags::ConstraintGamma1,
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, Dim, Frame>,
      ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<Dim, Frame>>,
      domain::Tags::MeshVelocity<Dim, Frame>>;

  using return_type = typename base::type;

  static void function(
      const gsl::not_null<return_type*> result,
      //   const Scalar<DataVector>& gamma_1,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, Dim, Frame>& shift,
      const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form,
      const std::optional<tnsr::I<DataVector, Dim, Frame>>& mesh_velocity) {
    characteristic_speeds(result,
                          // gamma_1,
                          lapse, shift, unit_normal_one_form, mesh_velocity);
  };
};

/// @}

/// @{
/*!
 * \brief Computes characteristic fields from evolved fields.
 */
template <size_t Dim, typename Frame>
typename Tags::CharacteristicFields<DataVector, Dim, Frame>::type
characteristic_fields(
    // const Scalar<DataVector>& gamma_2,
    const tnsr::II<DataVector, Dim, Frame>& inverse_spatial_metric,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    const tnsr::aa<DataVector, Dim, Frame>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame>& pi,
    // const tnsr::iaa<DataVector, Dim, Frame>& phi,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form);

template <size_t Dim, typename Frame>
void characteristic_fields(
    const gsl::not_null<
        typename Tags::CharacteristicFields<DataVector, Dim, Frame>::type*>
        char_fields,
    // const Scalar<DataVector>& gamma_2,
    const tnsr::II<DataVector, Dim, Frame>& inverse_spatial_metric,
    // Scalar driver fields
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    // Tensor driver fields
    const tnsr::aa<DataVector, Dim, Frame>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame>& pi,
    // const tnsr::iaa<DataVector, Dim, Frame>& phi,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form);

template <size_t Dim, typename Frame>
void characteristic_fields(
    const gsl::not_null<Scalar<DataVector>*>& v_scalar_driver,
    const gsl::not_null<Scalar<DataVector>*>& v_pi_scalar,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame>*>& v_tensor_driver,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame>*>& v_pi,
    // const Scalar<DataVector>& gamma_2,
    // const tnsr::II<DataVector, Dim, Frame>& inverse_spatial_metric,
    // Scalar driver fields
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    // Tensor driver fields
    const tnsr::aa<DataVector, Dim, Frame>& spacetime_metric,
    const tnsr::aa<DataVector, Dim, Frame>& pi,
    // const tnsr::iaa<DataVector, Dim, Frame>& phi,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form);

template <size_t Dim, typename Frame>
struct CharacteristicFieldsCompute
    : Tags::CharacteristicFields<DataVector, Dim, Frame>,
      db::ComputeTag {
  using base = Tags::CharacteristicFields<DataVector, Dim, Frame>;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<
      //   ::gh::ConstraintDamping::Tags::ConstraintGamma2,
      gr::Tags::InverseSpatialMetric<DataVector, Dim, Frame>, Tags::Psi,
      Tags::PiScalar, Tags::TensorDriver<DataVector, Dim, Frame>,
      Tags::Pi<DataVector, Dim, Frame>,
      //   Tags::Phi<DataVector, Dim, Frame>,
      ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<Dim, Frame>>>;

  static constexpr auto function = static_cast<void (*)(
      const gsl::not_null<return_type*>,
      const tnsr::II<DataVector, Dim, Frame>&, const Scalar<DataVector>&,
      const Scalar<DataVector>&, const tnsr::aa<DataVector, Dim, Frame>&,
      const tnsr::aa<DataVector, Dim, Frame>&,
      const tnsr::i<DataVector, Dim, Frame>&)>(
      &characteristic_fields<Dim, Frame>);
};
/// @}

/// @{
/*!
 * \brief For expressions used here to compute evolved fields from
 * characteristic ones, see \ref CharacteristicFieldsCompute.
 */
template <size_t Dim, typename Frame>
typename Tags::EvolvedFieldsFromCharacteristicFields<DataVector, Dim,
                                                     Frame>::type
evolved_fields_from_characteristic_fields(
    // const Scalar<DataVector>& gamma_2,
    const Scalar<DataVector>& u_scalar_driver,
    const Scalar<DataVector>& u_pi_scalar,
    const tnsr::aa<DataVector, Dim, Frame>& u_tensor_driver,
    const tnsr::aa<DataVector, Dim, Frame>& u_pi,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form);

template <size_t Dim, typename Frame>
void evolved_fields_from_characteristic_fields(
    const gsl::not_null<typename Tags::EvolvedFieldsFromCharacteristicFields<
        DataVector, Dim, Frame>::type*>
        evolved_fields,
    // const Scalar<DataVector>& gamma_2,
    const Scalar<DataVector>& u_scalar_driver,
    const Scalar<DataVector>& u_pi_scalar,
    const tnsr::aa<DataVector, Dim, Frame>& u_tensor_driver,
    const tnsr::aa<DataVector, Dim, Frame>& u_pi,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form);

template <size_t Dim, typename Frame>
struct EvolvedFieldsFromCharacteristicFieldsCompute
    : Tags::EvolvedFieldsFromCharacteristicFields<DataVector, Dim, Frame>,
      db::ComputeTag {
  using base =
      Tags::EvolvedFieldsFromCharacteristicFields<DataVector, Dim, Frame>;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<
      //   gh::ConstraintDamping::Tags::ConstraintGamma2,
      Tags::VScalarDriver<DataVector>, Tags::VPiScalar<DataVector>,
      Tags::VTensorDriver<DataVector, Dim, Frame>,
      Tags::VPi<DataVector, Dim, Frame>,
      ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<Dim, Frame>>>;

  static constexpr auto function = static_cast<void (*)(
      const gsl::not_null<return_type*>,
      const Scalar<DataVector>& u_scalar_driver,
      const Scalar<DataVector>& u_pi_scalar,
      const tnsr::aa<DataVector, Dim, Frame>& u_tensor_driver,
      const tnsr::aa<DataVector, Dim, Frame>& u_pi,
      const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form)>(
      &evolved_fields_from_characteristic_fields);
};
/// @}

namespace Tags {
struct LargestCharacteristicSpeed : db::SimpleTag {
  using type = double;
};
/*!
 * \brief Computes the largest magnitude of the characteristic speeds.
 */
// template <size_t Dim, typename Frame>
struct ComputeLargestCharacteristicSpeed : db::ComputeTag,
                                           LargestCharacteristicSpeed {
  using argument_tags = tmpl::list<
      //   ::gh::ConstraintDamping::Tags::ConstraintGamma1,
      gr::Tags::Lapse<DataVector>,
      gr::Tags::Shift<DataVector, 3, Frame::Inertial>,
      gr::Tags::SpatialMetric<DataVector, 3, Frame::Inertial>>;
  using return_type = double;
  using base = LargestCharacteristicSpeed;
  static void function(
      const gsl::not_null<double*> speed,
      //    const Scalar<DataVector>& gamma_1,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3, Frame::Inertial>& shift,
      const tnsr::ii<DataVector, 3, Frame::Inertial>& spatial_metric);
};
}  // namespace Tags
}  // namespace fe::ScalarTensorDriver
