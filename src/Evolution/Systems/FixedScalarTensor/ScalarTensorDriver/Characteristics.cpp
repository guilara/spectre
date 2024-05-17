// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Characteristics.hpp"

#include <algorithm>  // IWYU pragma: keep
#include <array>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"  // IWYU pragma: keep
#include "Domain/TagsTimeDependent.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

// IWYU pragma: no_forward_declare Tensor

namespace fe::ScalarTensorDriver {

template <size_t Dim, typename Frame>
void characteristic_speeds(
    const gsl::not_null<tnsr::a<DataVector, 3, Frame>*> char_speeds,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame>& shift,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form,
    const std::optional<tnsr::I<DataVector, Dim, Frame>>& mesh_velocity) {
  const auto shift_dot_normal = get(dot_product(shift, unit_normal_one_form));
  get<0>(*char_speeds) = -shift_dot_normal;  // lambda(VScalarDriver)
  get<1>(*char_speeds) = -shift_dot_normal;  // lambda(VPiScalar)
  get<2>(*char_speeds) = -shift_dot_normal;  // lambda(VTensorDriver)
  get<3>(*char_speeds) = -shift_dot_normal;  // lambda(VPi)
}

template <size_t Dim, typename Frame>
void characteristic_speeds(
    const gsl::not_null<std::array<DataVector, 4>*> char_speeds,
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame>& shift,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form,
    const std::optional<tnsr::I<DataVector, Dim, Frame>>& mesh_velocity) {
  const auto shift_dot_normal = get(dot_product(shift, unit_normal_one_form));
  (*char_speeds)[0] = -shift_dot_normal;  // lambda(VScalarDriver)
  (*char_speeds)[1] = -shift_dot_normal;  // lambda(VPiScalar)
  (*char_speeds)[2] = -shift_dot_normal;  // lambda(VTensorDriver)
  (*char_speeds)[3] = -shift_dot_normal;  // lambda(VPi)
}

template <size_t Dim, typename Frame>
std::array<DataVector, 4> characteristic_speeds(
    const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame>& shift,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form,
    const std::optional<tnsr::I<DataVector, Dim, Frame>>& mesh_velocity) {
  auto char_speeds = make_with_value<
      typename Tags::CharacteristicSpeeds<DataVector, Dim, Frame>::type>(
      get(lapse), 0.);
  characteristic_speeds(make_not_null(&char_speeds),

                        lapse, shift, unit_normal_one_form, mesh_velocity);
  return char_speeds;
}

template <size_t Dim, typename Frame>
void characteristic_fields(
    const gsl::not_null<
        typename Tags::CharacteristicFields<DataVector, Dim, Frame>::type*>
        char_fields,
    // Scalar driver fields
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    // Tensor driver fields
    const tnsr::aa<DataVector, Dim, Frame>& tensor_driver,
    const tnsr::aa<DataVector, Dim, Frame>& pi,

    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form) {
  const auto number_of_grid_points = get(psi).size();
  if (UNLIKELY(number_of_grid_points != char_fields->number_of_grid_points())) {
    char_fields->initialize(number_of_grid_points);
  }

  get<Tags::VScalarDriver<DataVector>>(*char_fields) = psi;
  get<Tags::VPiScalar<DataVector>>(*char_fields) = pi_scalar;
  get<Tags::VTensorDriver<DataVector, Dim, Frame>>(*char_fields) =
      tensor_driver;
  get<Tags::VPi<DataVector, Dim, Frame>>(*char_fields) = pi;
}

template <size_t Dim, typename Frame>
void characteristic_fields(
    const gsl::not_null<Scalar<DataVector>*>& v_scalar_driver,
    const gsl::not_null<Scalar<DataVector>*>& v_pi_scalar,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame>*>& v_tensor_driver,
    const gsl::not_null<tnsr::aa<DataVector, Dim, Frame>*>& v_pi,
    // Scalar driver fields
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    // Tensor driver fields
    const tnsr::aa<DataVector, Dim, Frame>& tensor_driver,
    const tnsr::aa<DataVector, Dim, Frame>& pi,

    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form) {
  *v_scalar_driver = psi;
  *v_pi_scalar = pi_scalar;
  *v_tensor_driver = tensor_driver;
  *v_pi = pi;
}

template <size_t Dim, typename Frame>
typename Tags::CharacteristicFields<DataVector, Dim, Frame>::type
characteristic_fields(
    const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
    const tnsr::aa<DataVector, Dim, Frame>& tensor_driver,
    const tnsr::aa<DataVector, Dim, Frame>& pi,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form) {
  auto char_fields = make_with_value<
      typename Tags::CharacteristicFields<DataVector, Dim, Frame>::type>(
      get(pi_scalar), 0.);
  characteristic_fields(make_not_null(&char_fields),

                        psi, pi_scalar, tensor_driver, pi,

                        unit_normal_one_form);
  return char_fields;
}

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
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form) {
  const auto number_of_grid_points = get(u_scalar_driver).size();
  if (UNLIKELY(number_of_grid_points !=
               evolved_fields->number_of_grid_points())) {
    evolved_fields->initialize(number_of_grid_points);
  }

  get<Tags::Psi>(*evolved_fields) = u_scalar_driver;
  get<Tags::PiScalar>(*evolved_fields) = u_pi_scalar;
  get<Tags::TensorDriver<DataVector, Dim, Frame>>(*evolved_fields) =
      u_tensor_driver;
  get<Tags::Pi<DataVector, Dim, Frame>>(*evolved_fields) = u_pi;
}

template <size_t Dim, typename Frame>
typename Tags::EvolvedFieldsFromCharacteristicFields<DataVector, Dim,
                                                     Frame>::type
evolved_fields_from_characteristic_fields(
    const Scalar<DataVector>& u_scalar_driver,
    const Scalar<DataVector>& u_pi_scalar,
    const tnsr::aa<DataVector, Dim, Frame>& u_tensor_driver,
    const tnsr::aa<DataVector, Dim, Frame>& u_pi,
    const tnsr::i<DataVector, Dim, Frame>& unit_normal_one_form) {
  auto evolved_fields =
      make_with_value<typename Tags::EvolvedFieldsFromCharacteristicFields<
          DataVector, Dim, Frame>::type>(get(u_scalar_driver), 0.);
  evolved_fields_from_characteristic_fields(make_not_null(&evolved_fields),
                                            // gamma_2,
                                            u_scalar_driver, u_pi_scalar,
                                            u_tensor_driver, u_pi,
                                            unit_normal_one_form);
  return evolved_fields;
}

template <size_t Dim, typename Frame>
void Tags::ComputeLargestCharacteristicSpeed<Dim, Frame>::function(
    const gsl::not_null<double*> speed, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, Dim, Frame>& shift,
    const tnsr::ii<DataVector, Dim, Frame>& spatial_metric) {
  const auto shift_magnitude = magnitude(shift, spatial_metric);
  *speed = max(get(shift_magnitude));
}
}  // namespace fe::ScalarTensorDriver

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATION(_, data)                                                 \
  template void fe::ScalarTensorDriver::characteristic_speeds(                 \
      const gsl::not_null<tnsr::a<DataVector, DIM(data), FRAME(data)>*>        \
          char_speeds,                                                         \
      const Scalar<DataVector>& lapse,                                         \
      const tnsr::I<DataVector, DIM(data), FRAME(data)>& shift,                \
      const tnsr::i<DataVector, DIM(data), FRAME(data)>& unit_normal_one_form, \
      const std::optional<tnsr::I<DataVector, DIM(data), FRAME(data)>>&        \
          mesh_velocity);                                                      \
  template void fe::ScalarTensorDriver::characteristic_speeds(                 \
      const gsl::not_null<std::array<DataVector, 4>*> char_speeds,             \
      const Scalar<DataVector>& lapse,                                         \
      const tnsr::I<DataVector, DIM(data), FRAME(data)>& shift,                \
      const tnsr::i<DataVector, DIM(data), FRAME(data)>& unit_normal_one_form, \
      const std::optional<tnsr::I<DataVector, DIM(data), FRAME(data)>>&        \
          mesh_velocity);                                                      \
  template std::array<DataVector, 4>                                           \
  fe::ScalarTensorDriver::characteristic_speeds(                               \
      const Scalar<DataVector>& lapse,                                         \
      const tnsr::I<DataVector, DIM(data), FRAME(data)>& shift,                \
      const tnsr::i<DataVector, DIM(data), FRAME(data)>& unit_normal_one_form, \
      const std::optional<tnsr::I<DataVector, DIM(data), FRAME(data)>>&        \
          mesh_velocity);                                                      \
  template struct fe::ScalarTensorDriver::CharacteristicSpeedsCompute<         \
      DIM(data), FRAME(data)>;                                                 \
  template void fe::ScalarTensorDriver::characteristic_fields(                 \
      const gsl::not_null<                                                     \
          typename fe::ScalarTensorDriver::Tags::CharacteristicFields<         \
              DataVector, DIM(data), FRAME(data)>::type*>                      \
          char_fields,                                                         \
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,      \
      const tnsr::aa<DataVector, DIM(data), FRAME(data)>& tensor_driver,       \
      const tnsr::aa<DataVector, DIM(data), FRAME(data)>& pi,                  \
      const tnsr::i<DataVector, DIM(data), FRAME(data)>&                       \
          unit_normal_one_form);                                               \
  template void fe::ScalarTensorDriver::characteristic_fields(                 \
      const gsl::not_null<Scalar<DataVector>*>& v_scalar_driver,               \
      const gsl::not_null<Scalar<DataVector>*>& v_pi_scalar,                   \
      const gsl::not_null<tnsr::aa<DataVector, DIM(data), FRAME(data)>*>&      \
          v_tensor_driver,                                                     \
      const gsl::not_null<tnsr::aa<DataVector, DIM(data), FRAME(data)>*>&      \
          v_pi,                                                                \
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,      \
      const tnsr::aa<DataVector, DIM(data), FRAME(data)>& tensor_driver,       \
      const tnsr::aa<DataVector, DIM(data), FRAME(data)>& pi,                  \
      const tnsr::i<DataVector, DIM(data), FRAME(data)>&                       \
          unit_normal_one_form);                                               \
  template typename fe::ScalarTensorDriver::Tags::CharacteristicFields<        \
      DataVector, DIM(data), FRAME(data)>::type                                \
  fe::ScalarTensorDriver::characteristic_fields(                               \
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,      \
      const tnsr::aa<DataVector, DIM(data), FRAME(data)>& tensor_driver,       \
      const tnsr::aa<DataVector, DIM(data), FRAME(data)>& pi,                  \
      const tnsr::i<DataVector, DIM(data), FRAME(data)>&                       \
          unit_normal_one_form);                                               \
  template struct fe::ScalarTensorDriver::CharacteristicFieldsCompute<         \
      DIM(data), FRAME(data)>;                                                 \
  template void                                                                \
  fe::ScalarTensorDriver::evolved_fields_from_characteristic_fields(           \
      const gsl::not_null<typename fe::ScalarTensorDriver::Tags::              \
                              EvolvedFieldsFromCharacteristicFields<           \
                                  DataVector, DIM(data), FRAME(data)>::type*>  \
          evolved_fields,                                                      \
      const Scalar<DataVector>& u_scalar_driver,                               \
      const Scalar<DataVector>& u_pi_scalar,                                   \
      const tnsr::aa<DataVector, DIM(data), FRAME(data)>& u_tensor_driver,     \
      const tnsr::aa<DataVector, DIM(data), FRAME(data)>& u_pi,                \
      const tnsr::i<DataVector, DIM(data), FRAME(data)>&                       \
          unit_normal_one_form);                                               \
  template typename fe::ScalarTensorDriver::Tags::                             \
      EvolvedFieldsFromCharacteristicFields<DataVector, DIM(data),             \
                                            FRAME(data)>::type                 \
      fe::ScalarTensorDriver::evolved_fields_from_characteristic_fields(       \
          const Scalar<DataVector>& u_scalar_driver,                           \
          const Scalar<DataVector>& u_pi_scalar,                               \
          const tnsr::aa<DataVector, DIM(data), FRAME(data)>& u_tensor_driver, \
          const tnsr::aa<DataVector, DIM(data), FRAME(data)>& u_pi,            \
          const tnsr::i<DataVector, DIM(data), FRAME(data)>&                   \
              unit_normal_one_form);                                           \
  template struct fe::ScalarTensorDriver::                                     \
      EvolvedFieldsFromCharacteristicFieldsCompute<DIM(data), FRAME(data)>;    \
  template struct fe::ScalarTensorDriver::Tags::                               \
      ComputeLargestCharacteristicSpeed<DIM(data), FRAME(data)>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (3),
                        (Frame::Inertial, Frame::Grid))

#undef INSTANTIATION
#undef DIM
#undef FRAME
