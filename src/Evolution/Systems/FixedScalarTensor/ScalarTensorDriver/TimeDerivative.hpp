// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedContainers.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "Domain/TagsTimeDependent.hpp"
#include "Evolution/Systems/CurvedScalarWave/TimeDerivative.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::ScalarTensorDriver {
/*!
 * \brief Compute the RHS terms of the evolution equations for the scalar tensor
 * driver system.
 */
struct TimeDerivative {
  static constexpr size_t dim = 3;

  using dt_tags =
      db::wrap_tags_in<::Tags::dt, typename System::variables_tag::tags_list>;

  using gradient_tags = typename System::variables_tag::gradients_tags;

  using temporary_tags = tmpl::list<
      // Tensor driver temporary tags
      ::gh::ConstraintDamping::Tags::ConstraintGamma1,
      ::gh::ConstraintDamping::Tags::ConstraintGamma2,
      // gh::Tags::GaugeH<DataVector, Dim>,
      // gh::Tags::SpacetimeDerivGaugeH<DataVector, Dim>,
      // gh::Tags::Gamma1Gamma2, gh::Tags::HalfPiTwoNormals,
      // gh::Tags::NormalDotOneIndexConstraint, gh::Tags::Gamma1Plus1,
      // gh::Tags::PiOneNormal<Dim>, gh::Tags::GaugeConstraint<DataVector, Dim>,
      // gh::Tags::HalfPhiTwoNormals<Dim>,
      // gh::Tags::ShiftDotThreeIndexConstraint<Dim>,
      // gh::Tags::MeshVelocityDotThreeIndexConstraint<Dim>,
      // gh::Tags::PhiOneNormal<Dim>, gh::Tags::PiSecondIndexUp<Dim>,
      // gh::Tags::ThreeIndexConstraint<DataVector, Dim>,
      // gh::Tags::PhiFirstIndexUp<Dim>, gh::Tags::PhiThirdIndexUp<Dim>,
      // gh::Tags::SpacetimeChristoffelFirstKindThirdIndexUp<Dim>,
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, Dim>,
      gr::Tags::InverseSpatialMetric<DataVector, Dim>,
      gr::Tags::DetSpatialMetric<DataVector>,
      // gr::Tags::SqrtDetSpatialMetric<DataVector>,
      // gr::Tags::InverseSpacetimeMetric<DataVector, Dim>,
      // gr::Tags::SpacetimeChristoffelFirstKind<DataVector, Dim>,
      // gr::Tags::SpacetimeChristoffelSecondKind<DataVector, Dim>,
      // gr::Tags::TraceSpacetimeChristoffelFirstKind<DataVector, Dim>,
      // gr::Tags::SpacetimeNormalVector<DataVector, Dim>,
      // Scalar driver temporary tags
      // gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, Dim>,
      // gr::Tags::InverseSpatialMetric<DataVector, Dim>,
      CurvedScalarWave::Tags::ConstraintGamma1,
      CurvedScalarWave::Tags::ConstraintGamma2,
      // Extra scalar driver temporary tags
      ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, dim,
                                                    ::Frame::Inertial>>;

  using argument_tags = tmpl::list<
      // Tensor argument tags
      gr::Tags::SpacetimeMetric<DataVector, Dim>, gh::Tags::Pi<DataVector, Dim>,
      gh::Tags::Phi<DataVector, Dim>,
      ::gh::ConstraintDamping::Tags::ConstraintGamma0,
      ::gh::ConstraintDamping::Tags::ConstraintGamma1,
      ::gh::ConstraintDamping::Tags::ConstraintGamma2,
      gauges::Tags::GaugeCondition, domain::Tags::Mesh<Dim>, ::Tags::Time,
      domain::Tags::Coordinates<Dim, Frame::Inertial>,
      domain::Tags::InverseJacobian<Dim, Frame::ElementLogical,
                                    Frame::Inertial>,
      domain::Tags::MeshVelocity<Dim, Frame::Inertial>,
      // Scalar argument tags
      CurvedScalarWave::Tags::Pi, CurvedScalarWave::Tags::Phi<Dim>,
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, Dim>,
      ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<Dim>,
                    Frame::Inertial>,
      ::Tags::deriv<gr::Tags::Shift<DataVector, Dim>, tmpl::size_t<Dim>,
                    Frame::Inertial>,
      gr::Tags::InverseSpatialMetric<DataVector, Dim>,
      gr::Tags::TraceSpatialChristoffelSecondKind<DataVector, Dim>,
      gr::Tags::TraceExtrinsicCurvature<DataVector>,
      CurvedScalarWave::Tags::ConstraintGamma1,
      CurvedScalarWave::Tags::ConstraintGamma2,
      // Extra argument tags
      ScalarTensor::Tags::ScalarSource>;

  static void apply(
      // GH dt variables
      gsl::not_null<tnsr::aa<DataVector, dim>*> dt_spacetime_metric,
      gsl::not_null<tnsr::aa<DataVector, dim>*> dt_pi,
      gsl::not_null<tnsr::iaa<DataVector, dim>*> dt_phi,
      // Scalar dt variables
      gsl::not_null<Scalar<DataVector>*> dt_psi_scalar,
      gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
      gsl::not_null<tnsr::i<DataVector, dim, Frame::Inertial>*> dt_phi_scalar,

      // GH temporal variables
      gsl::not_null<Scalar<DataVector>*> temp_gamma1,
      gsl::not_null<Scalar<DataVector>*> temp_gamma2,
      // gsl::not_null<tnsr::a<DataVector, dim>*> temp_gauge_function,
      // gsl::not_null<tnsr::ab<DataVector, dim>*>
      //     temp_spacetime_deriv_gauge_function,
      // gsl::not_null<Scalar<DataVector>*> gamma1gamma2,
      // gsl::not_null<Scalar<DataVector>*> half_half_pi_two_normals,
      // gsl::not_null<Scalar<DataVector>*> normal_dot_gauge_constraint,
      // gsl::not_null<Scalar<DataVector>*> gamma1_plus_1,
      // gsl::not_null<tnsr::a<DataVector, dim>*> pi_one_normal,
      // gsl::not_null<tnsr::a<DataVector, dim>*> gauge_constraint,
      // gsl::not_null<tnsr::i<DataVector, dim>*> half_phi_two_normals,
      // gsl::not_null<tnsr::aa<DataVector, dim>*>
      //     shift_dot_three_index_constraint,
      // gsl::not_null<tnsr::aa<DataVector, dim>*>
      //     mesh_velocity_dot_three_index_constraint,
      // gsl::not_null<tnsr::ia<DataVector, dim>*> phi_one_normal,
      // gsl::not_null<tnsr::aB<DataVector, dim>*> pi_2_up,
      // gsl::not_null<tnsr::iaa<DataVector, dim>*> three_index_constraint,
      // gsl::not_null<tnsr::Iaa<DataVector, dim>*> phi_1_up,
      // gsl::not_null<tnsr::iaB<DataVector, dim>*> phi_3_up,
      // gsl::not_null<tnsr::abC<DataVector, dim>*> christoffel_first_kind_3_up,
      gsl::not_null<Scalar<DataVector>*> lapse,
      gsl::not_null<tnsr::I<DataVector, dim>*> shift,
      gsl::not_null<tnsr::II<DataVector, dim>*> inverse_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> det_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> sqrt_det_spatial_metric,
      // gsl::not_null<tnsr::AA<DataVector, dim>*> inverse_spacetime_metric,
      // gsl::not_null<tnsr::abb<DataVector, dim>*> christoffel_first_kind,
      // gsl::not_null<tnsr::Abb<DataVector, dim>*> christoffel_second_kind,
      // gsl::not_null<tnsr::a<DataVector, dim>*> trace_christoffel,
      // gsl::not_null<tnsr::A<DataVector, dim>*> normal_spacetime_vector,

      // Scalar temporal variables
      gsl::not_null<Scalar<DataVector>*> result_gamma1_scalar,
      gsl::not_null<Scalar<DataVector>*> result_gamma2_scalar,

      // Extra temporal tags
      gsl::not_null<tnsr::aa<DataVector, dim>*> stress_energy,

      // GH spatial derivatives
      const tnsr::iaa<DataVector, dim>& d_spacetime_metric,
      const tnsr::iaa<DataVector, dim>& d_pi,
      const tnsr::ijaa<DataVector, dim>& d_phi,

      // scalar spatial derivatives
      const tnsr::i<DataVector, dim>& d_psi_scalar,
      const tnsr::i<DataVector, dim>& d_pi_scalar,
      const tnsr::ij<DataVector, dim>& d_phi_scalar,

      // GH argument variables
      const tnsr::aa<DataVector, dim>& spacetime_metric,
      const tnsr::aa<DataVector, dim>& pi,
      const tnsr::iaa<DataVector, dim>& phi, const Scalar<DataVector>& gamma0,
      const Scalar<DataVector>& gamma1, const Scalar<DataVector>& gamma2,
      const gh::gauges::GaugeCondition& gauge_condition, const Mesh<dim>& mesh,
      double time,
      const tnsr::I<DataVector, dim, Frame::Inertial>& inertial_coords,
      const InverseJacobian<DataVector, dim, Frame::ElementLogical,
                            Frame::Inertial>& inverse_jacobian,
      const std::optional<tnsr::I<DataVector, dim, Frame::Inertial>>&
          mesh_velocity,

      // Scalar argument variables
      const Scalar<DataVector>& pi_scalar,
      const tnsr::i<DataVector, dim>& phi_scalar,
      const Scalar<DataVector>& lapse_scalar,
      const tnsr::I<DataVector, dim>& shift_scalar,
      const tnsr::i<DataVector, dim>& deriv_lapse,
      const tnsr::iJ<DataVector, dim>& deriv_shift,
      const tnsr::II<DataVector, dim>& upper_spatial_metric,
      const tnsr::I<DataVector, dim>& trace_spatial_christoffel,
      const Scalar<DataVector>& trace_extrinsic_curvature,
      const Scalar<DataVector>& gamma1_scalar,
      const Scalar<DataVector>& gamma2_scalar,

      const Scalar<DataVector>& scalar_source);
};
}  // namespace fe::ScalarTensorDriver
