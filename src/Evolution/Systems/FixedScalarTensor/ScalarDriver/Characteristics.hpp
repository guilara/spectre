// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/FaceNormal.hpp"
#include "Evolution/Systems/CurvedScalarWave/Characteristics.hpp"
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

  static void function(
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
    const tnsr::I<DataVector, 3_st, Frame::Inertial>& unit_normal_vector);

struct CharacteristicFieldsCompute : Tags::CharacteristicFields,
                                     db::ComputeTag {
  using base = Tags::CharacteristicFields;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<
      Tags::ConstraintGamma2, gr::Tags::InverseSpatialMetric<DataVector, 3_st>,
      Tags::Psi, Tags::Pi, Tags::Phi<3_st>,
      ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<3_st>>>;

  static void function(
      gsl::not_null<return_type*> result, const Scalar<DataVector>& gamma_2,
      const tnsr::II<DataVector, 3_st, Frame::Inertial>& inverse_spatial_metric,
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& phi,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form) {
    const auto unit_normal_vector = tenex::evaluate<ti::I>(
        inverse_spatial_metric(ti::I, ti::J) * unit_normal_one_form(ti::j));
    characteristic_fields(result, gamma_2, psi, pi, phi, unit_normal_one_form,
                          unit_normal_vector);
  }
};

Variables<tmpl::list<Tags::Psi, Tags::Pi, Tags::Phi<3_st>>>
evolved_fields_from_characteristic_fields(
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& v_psi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero,
    const Scalar<DataVector>& v_plus, const Scalar<DataVector>& v_minus,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form);

void evolved_fields_from_characteristic_fields(
    gsl::not_null<Variables<tmpl::list<Tags::Psi, Tags::Pi, Tags::Phi<3_st>>>*>
        evolved_fields,
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& v_psi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero,
    const Scalar<DataVector>& v_plus, const Scalar<DataVector>& v_minus,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form);

void evolved_fields_from_characteristic_fields(
    gsl::not_null<Scalar<DataVector>*> psi,
    gsl::not_null<Scalar<DataVector>*> pi,
    gsl::not_null<tnsr::i<DataVector, 3_st, Frame::Inertial>*> phi,
    const Scalar<DataVector>& gamma_2, const Scalar<DataVector>& v_psi,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero,
    const Scalar<DataVector>& v_plus, const Scalar<DataVector>& v_minus,
    const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form);

struct EvolvedFieldsFromCharacteristicFieldsCompute
    : Tags::EvolvedFieldsFromCharacteristicFields,
      db::ComputeTag {
  using base = Tags::EvolvedFieldsFromCharacteristicFields;
  using return_type = typename base::type;
  using argument_tags = tmpl::list<
      Tags::ConstraintGamma2, Tags::VPsi, Tags::VZero<3_st>, Tags::VPlus,
      Tags::VMinus,
      ::Tags::Normalized<domain::Tags::UnnormalizedFaceNormal<3_st>>>;

  static void function(
      gsl::not_null<return_type*> result, const Scalar<DataVector>& gamma_2,
      const Scalar<DataVector>& v_psi,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& v_zero,
      const Scalar<DataVector>& v_plus, const Scalar<DataVector>& v_minus,
      const tnsr::i<DataVector, 3_st, Frame::Inertial>& unit_normal_one_form) {
    evolved_fields_from_characteristic_fields(
        result, gamma_2, v_psi, v_zero, v_plus, v_minus, unit_normal_one_form);
  }
};

namespace Tags {

struct ComputeLargestCharacteristicSpeed : LargestCharacteristicSpeed,
                                           db::ComputeTag {
  using argument_tags =
      tmpl::list<Tags::ConstraintGamma1, gr::Tags::Lapse<DataVector>,
                 gr::Tags::Shift<DataVector, 3_st>,
                 gr::Tags::SpatialMetric<DataVector, 3_st>>;
  using return_type = double;
  using base = LargestCharacteristicSpeed;
  static void function(
      const gsl::not_null<double*> max_speed, const Scalar<DataVector>& gamma_1,
      const Scalar<DataVector>& lapse,
      const tnsr::I<DataVector, 3_st, Frame::Inertial>& shift,
      const tnsr::ii<DataVector, 3_st, Frame::Inertial>& spatial_metric) {
    // Use the methods of CurvedScalarWave
    CurvedScalarWave::Tags::ComputeLargestCharacteristicSpeed<3_st>::function(
        max_speed, gamma_1, lapse, shift, spatial_metric);
  }
};

}  // namespace Tags

}  // namespace fe::ScalarDriver
