// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/FaceNormal.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/TagsDeclarations.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::ScalarDriver {

struct CharacteristicSpeedsCompute : Tags::CharacteristicSpeeds,
                                     db::ComputeTag {
  using base = Tags::CharacteristicSpeeds;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<
      Tags::ConstraintGamma1, gr::Tags::Lapse<DataVector>,
      gr::Tags::Shift<DataVector, 3_st>,
      ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<3_st>>>;

  static constexpr void function(
      gsl::not_null<return_type*> result, const Scalar<DataVector>& gamma_1,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form) {
    CurvedScalarWave::characteristic_speeds<3_st>(result, gamma_1, lapse, shift,
                                                  unit_normal_one_form);
  }
};

Variables<tmpl::list<Tags::VPsi, Tags::VZero<3_st>, Tags::VPlus, Tags::VMinus>>
characteristic_fields(
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& psi,
    const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>& unit_normal_vector);

void characteristic_fields(
    gsl::not_null<Variables<
        tmpl::list<Tags::VPsi, Tags::VZero<3_st>, Tags::VPlus, Tags::VMinus>>*>
        char_fields,
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& psi,
    const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>& unit_normal_vector);

void characteristic_fields(
    const gsl::not_null<Scalar<DataVector>*>& v_psi,
    const gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*>& v_zero,
    const gsl::not_null<Scalar<DataVector>*>& v_plus,
    const gsl::not_null<Scalar<DataVector>*>& v_minus,
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& psi,
    const Scalar<DataVector>& pi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form,
    const tnsr::I<DataVector, 3_st, Frame::Inertial>& unit_normal_vector) {
  // Use the methods from CurvedScalarWave
  CurvedScalarWave::characteristic_fields<3_st>(
      v_psi, v_zero, v_plus, v_minus, psi, pi, phi, unit_normal_one_form,
      unit_normal_vector);
}

namespace Tags {}  // namespace Tags

}  // namespace fe::ScalarDriver
