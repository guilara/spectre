// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Characteristics.hpp"

#include <algorithm>
#include <array>
#include <limits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace fe::ScalarDriver {

void characteristic_speeds(
    const gsl::not_null<tnsr::a<DataVector, 3, Frame::Inertial>*> char_speeds,
    const Scalar<DataVector>& gamma_1, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form) {
  destructive_resize_components(char_speeds, get(gamma_1).size());
  const auto shift_dot_normal = get(dot_product(shift, unit_normal_one_form));
  get<0>(*char_speeds) = -(1. + get(gamma_1)) * shift_dot_normal;  // v(VPsi)
  get<1>(*char_speeds) = -shift_dot_normal;                        // v(VZero)
  get<2>(*char_speeds) = -shift_dot_normal;                        // v(VPlus)
  get<3>(*char_speeds) = -shift_dot_normal;                        // v(VMinus)
}

void characteristic_speeds(
    const gsl::not_null<std::array<DataVector, 4>*> char_speeds,
    const Scalar<DataVector>& gamma_1, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form) {
  const size_t size = get(gamma_1).size();
  for (auto& char_speed : *char_speeds) {
    char_speed.destructive_resize(size);
  }
  const auto shift_dot_normal = get(dot_product(shift, unit_normal_one_form));
  (*char_speeds)[0] = -(1. + get(gamma_1)) * shift_dot_normal;  // v(VPsi)
  (*char_speeds)[1] = -shift_dot_normal;                        // v(VZero)
  (*char_speeds)[2] = -shift_dot_normal;                        // v(VPlus)
  (*char_speeds)[3] = -shift_dot_normal;                        // v(VMinus)
}

std::array<DataVector, 4> characteristic_speeds(
    const Scalar<DataVector>& gamma_1, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form) {
  std::array<DataVector, 4> char_speeds{};
  characteristic_speeds(make_not_null(&char_speeds), gamma_1, lapse, shift,
                        unit_normal_one_form);
  return char_speeds;
}

void characteristic_fields(
    const gsl::not_null<Variables<tmpl::list<
        Tags::VPsi, Tags::VZero<3_st>, Tags::VPlus, Tags::VMinus>>*>
        char_fields,
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& psi,
    const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>&
        unit_normal_one_form,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>&
        unit_normal_vector) {
  char_fields->initialize(get(gamma_2).size());
  dot_product(make_not_null(&get<Tags::VMinus>(*char_fields)),
              unit_normal_vector, phi);
  // Eq.(34) of Holst+ (2004) for VZero
  for (size_t i = 0; i < 3_st; ++i) {
    get<Tags::VZero<3_st>>(*char_fields).get(i) = 0.0 * phi.get(i);
  }
  // Eq.(33) of Holst+ (2004) for VPsi
  get(get<Tags::VPsi>(*char_fields)) = get(psi);
  // Eq.(35) of Holst+ (2004) for VPlus and VMinus
  get(get<Tags::VPlus>(*char_fields)) = get(pi);
  get(get<Tags::VMinus>(*char_fields)) =
      0.0 * get(get<Tags::VMinus>(*char_fields));
}

Variables<
    tmpl::list<Tags::VPsi, Tags::VZero<3_st>, Tags::VPlus, Tags::VMinus>>
characteristic_fields(
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& psi,
    const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>&
        unit_normal_one_form,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>&
        unit_normal_vector) {
  Variables<tmpl::list<Tags::VPsi, Tags::VZero<3_st>, Tags::VPlus,
                       Tags::VMinus>>
      char_fields{get(gamma_2).size()};
  characteristic_fields(make_not_null(&char_fields), gamma_2, psi, pi, phi,
                        unit_normal_one_form, unit_normal_vector);
  return char_fields;
}

void characteristic_fields(
    const gsl::not_null<Scalar<DataVector>*>& v_psi,
    const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>&
        v_zero,
    const gsl::not_null<Scalar<DataVector>*>& v_plus,
    const gsl::not_null<Scalar<DataVector>*>& v_minus,
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& psi,
    const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>&
        unit_normal_one_form,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>&
        unit_normal_vector) {
  const size_t size = get(gamma_2).size();
  destructive_resize_components(v_psi, size);
  destructive_resize_components(v_zero, size);
  destructive_resize_components(v_plus, size);
  destructive_resize_components(v_minus, size);

  dot_product(v_minus, unit_normal_vector, phi);
  // Eq.(34) of Holst+ (2004) for VZero
  for (size_t i = 0; i < 3_st; ++i) {
    v_zero->get(i) = 0.0 * phi.get(i);
  }
  // Eq.(33) of Holst+ (2004) for VPsi
  get(*v_psi) = get(psi);
  // Eq.(35) of Holst+ (2004) for VPlus and VMinus
  get(*v_plus) = get(pi);
  get(*v_minus) = 0.0 * get(*v_minus);
}

void evolved_fields_from_characteristic_fields(
    gsl::not_null<Scalar<DataVector>*> psi,
    gsl::not_null<Scalar<DataVector>*> pi,
    gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*> phi,
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& v_psi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero,
    const Scalar<DataVector>& v_plus, const Scalar<DataVector>& v_minus,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>&
        unit_normal_one_form) {
  const size_t size = get(gamma_2).size();
  destructive_resize_components(psi, size);
  destructive_resize_components(pi, size);
  destructive_resize_components(phi, size);
  // Eq.(36) of Holst+ (2005) for Psi
  psi->get() = get(v_psi);
  // Eq.(37) - (38) of Holst+ (2004) for Pi and Phi
  pi->get() = get(v_plus);
  for (size_t i = 0; i < 3_st; ++i) {
    phi->get(i) = get(v_minus) * unit_normal_one_form.get(i);
  }
}

void evolved_fields_from_characteristic_fields(
    const gsl::not_null<
        Variables<tmpl::list<Tags::Psi, Tags::Pi, Tags::Phi<3_st>>>*>
        evolved_fields,
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& v_psi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero,
    const Scalar<DataVector>& v_plus, const Scalar<DataVector>& v_minus,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>&
        unit_normal_one_form) {
  evolved_fields->initialize(get_size(get(gamma_2)));
  // Eq.(36) of Holst+ (2005) for Psi
  get<Tags::Psi>(*evolved_fields).get() = get(v_psi);

  // Eq.(37) - (38) of Holst+ (2004) for Pi and Phi
  get<Tags::Pi>(*evolved_fields).get() = get(v_plus);
  for (size_t i = 0; i < 3_st; ++i) {
    get<Tags::Phi<3_st>>(*evolved_fields).get(i) =
        get(v_minus) * unit_normal_one_form.get(i);
  }
}

Variables<tmpl::list<Tags::Psi, Tags::Pi, Tags::Phi<3_st>>>
evolved_fields_from_characteristic_fields(
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& v_psi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero,
    const Scalar<DataVector>& v_plus, const Scalar<DataVector>& v_minus,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>&
        unit_normal_one_form) {
  Variables<tmpl::list<Tags::Psi, Tags::Pi, Tags::Phi<3_st>>>
      evolved_fields(get(gamma_2).size());
  evolved_fields_from_characteristic_fields(make_not_null(&evolved_fields),
                                            gamma_2, v_psi, v_zero, v_plus,
                                            v_minus, unit_normal_one_form);
  return evolved_fields;
}

} // namespace fe::ScalarDriver
