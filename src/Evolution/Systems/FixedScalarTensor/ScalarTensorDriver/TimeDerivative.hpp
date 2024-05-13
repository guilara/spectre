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
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/System.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/InverseSpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Lapse.hpp"
#include "PointwiseFunctions/GeneralRelativity/Shift.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::ScalarTensorDriver {
// struct System;
/*!
 * \brief Compute the RHS terms of the evolution equations for the scalar tensor
 * driver system.
 */
struct TimeDerivative {
  static constexpr size_t dim = 3;

  using dt_tags = db::wrap_tags_in<
      ::Tags::dt,
      typename ::fe::ScalarTensorDriver::System::variables_tag::tags_list>;

  using gradient_tags =
      typename ::fe::ScalarTensorDriver::System::gradients_tags;

  using temporary_tags = tmpl::list<
      // Tensor driver temporary tags
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, dim>,
      gr::Tags::InverseSpatialMetric<DataVector, dim>,
      gr::Tags::DetSpatialMetric<DataVector>
      // Extra scalar driver temporary tags
      >;

  using argument_tags = tmpl::list<
      // Tensor argument tags
      fe::ScalarTensorDriver::Tags::TensorDriver<DataVector, dim>,
      fe::ScalarTensorDriver::Tags::Pi<DataVector, dim>,

      domain::Tags::Mesh<dim>, ::Tags::Time,
      domain::Tags::Coordinates<dim, Frame::Inertial>,
      domain::Tags::InverseJacobian<dim, Frame::ElementLogical,
                                    Frame::Inertial>,
      domain::Tags::MeshVelocity<dim, Frame::Inertial>,
      // Scalar argument tags
      fe::ScalarTensorDriver::Tags::Psi, fe::ScalarTensorDriver::Tags::PiScalar,
      gr::Tags::Lapse<DataVector>, gr::Tags::Shift<DataVector, dim>,

      // Extra argument tags
      fe::ScalarTensorDriver::Tags::TensorDriverSource<DataVector, dim,
                                                       ::Frame::Inertial>,
      fe::ScalarTensorDriver::Tags::ScalarDriverSource,
      // TODO: Add sigma and tau parameters
      fe::ScalarTensorDriver::Tags::TauParameter,
      fe::ScalarTensorDriver::Tags::SigmaParameter>;

  static void apply(
      // GH dt variables
      gsl::not_null<tnsr::aa<DataVector, dim>*> dt_spacetime_metric,
      gsl::not_null<tnsr::aa<DataVector, dim>*> dt_pi,
      //   gsl::not_null<tnsr::iaa<DataVector, dim>*> dt_phi,
      // Scalar dt variables
      gsl::not_null<Scalar<DataVector>*> dt_psi_scalar,
      gsl::not_null<Scalar<DataVector>*> dt_pi_scalar,
      //   gsl::not_null<tnsr::i<DataVector, dim, Frame::Inertial>*>
      //   dt_phi_scalar,

      // GH temporal variables
      gsl::not_null<Scalar<DataVector>*> lapse,
      gsl::not_null<tnsr::I<DataVector, dim>*> shift,
      gsl::not_null<tnsr::II<DataVector, dim>*> inverse_spatial_metric,
      gsl::not_null<Scalar<DataVector>*> det_spatial_metric,

      // Scalar temporal variables

      // Extra temporal tags
      //   gsl::not_null<tnsr::aa<DataVector, dim>*> stress_energy,

      // GH spatial derivatives
      const tnsr::iaa<DataVector, dim>& d_spacetime_metric,
      const tnsr::iaa<DataVector, dim>& d_pi,
      //   const tnsr::ijaa<DataVector, dim>& d_phi,

      // scalar spatial derivatives
      const tnsr::i<DataVector, dim>& d_psi_scalar,
      const tnsr::i<DataVector, dim>& d_pi_scalar,
      //   const tnsr::ij<DataVector, dim>& d_phi_scalar,

      // GH argument variables
      const tnsr::aa<DataVector, dim>& spacetime_metric,
      const tnsr::aa<DataVector, dim>& pi,

      const Mesh<dim>& mesh, double time,
      const tnsr::I<DataVector, dim, Frame::Inertial>& inertial_coords,
      const InverseJacobian<DataVector, dim, Frame::ElementLogical,
                            Frame::Inertial>& inverse_jacobian,
      const std::optional<tnsr::I<DataVector, dim, Frame::Inertial>>&
          mesh_velocity,

      // Scalar argument variables
      const Scalar<DataVector>& psi, const Scalar<DataVector>& pi_scalar,
      //   const tnsr::i<DataVector, dim>& phi_scalar,
      const Scalar<DataVector>& lapse_scalar,
      const tnsr::I<DataVector, dim>& shift_scalar,

      const tnsr::aa<DataVector, dim>& tensor_driver_source,
      const Scalar<DataVector>& scalar_driver_source,

      const Scalar<DataVector>& tau_parameter,
      const Scalar<DataVector>& sigma_parameter);
};
}  // namespace fe::ScalarTensorDriver
