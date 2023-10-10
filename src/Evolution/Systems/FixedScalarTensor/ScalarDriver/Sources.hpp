// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
//
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"
//
#include "Options/String.hpp"

namespace fe::ScalarDriver::Sources {

void add_scalar_driver_source_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi,
    const Scalar<DataVector>& scalar_driver_source,
    const Scalar<DataVector>& lapse);

void add_scalar_driver_friction_term_to_dt_pi_scalar(
    gsl::not_null<Scalar<DataVector>*> dt_pi,
    const Scalar<DataVector>& scalar_driver_pi, const Scalar<DataVector>& lapse,
    const tnsr::I<DataVector, 3_st>& shift, const double scalar_tau_parameter,
    const double scalar_sigma_parameter);

void compute_scalar_driver_source(gsl::not_null<Scalar<DataVector>*> result,
                                  const Scalar<DataVector>& psi,
                                  const Scalar<DataVector>& target_psi,
                                  const double scalar_tau_parameter,
                                  const double scalar_sigma_parameter);

void compute_scalar_driver_source_with_limiter(
    const gsl::not_null<Scalar<DataVector>*> scalar_driver_source,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const double scalar_tau_parameter, const double scalar_sigma_parameter,
    const double limiter_parameter);

void compute_target_psi(gsl::not_null<Scalar<DataVector>*> target_psi,
                        const Scalar<DataVector>& psi);

// For the exponential driver
void add_scalar_driver_source_to_dt_psi_scalar_for_exponential_driver(
    gsl::not_null<Scalar<DataVector>*> dt_psi,
    const Scalar<DataVector>& scalar_driver_source);

void compute_scalar_driver_source_for_exponential_driver(
    const gsl::not_null<Scalar<DataVector>*> scalar_driver_source,
    const Scalar<DataVector>& psi, const Scalar<DataVector>& target_psi,
    const double scalar_tau_parameter, const double scalar_sigma_parameter);

}  // namespace fe::ScalarDriver::Sources

namespace fe::ScalarDriver::Tags {

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <typename Frame, typename DataType>
struct ScalarDriverSourceCompute : ScalarDriverSource, db::ComputeTag {
  using argument_tags =
      tmpl::list<fe::ScalarDriver::Tags::Psi, fe::ScalarDriver::Tags::TargetPsi,
                 fe::ScalarDriver::Tags::ScalarTauParameter,
                 fe::ScalarDriver::Tags::ScalarSigmaParameter>;
  using return_type = Scalar<DataType>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataType>&,
      const Scalar<DataType>&, const double,
      const double) = &fe::ScalarDriver::Sources::
                          compute_scalar_driver_source_for_exponential_driver;
  using base = ScalarDriverSource;
};

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <typename Frame, typename DataType>
struct TargetPsiCompute : TargetPsi, db::ComputeTag {
  //   using argument_tags = tmpl::list<fe::ScalarDriver::Psi>;
  using argument_tags =
      tmpl::list<gr::Tags::WeylElectricScalar<DataType>,
                 gr::Tags::WeylMagneticScalar<DataType>,
                 CurvedScalarWave::Tags::Psi,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter,
                 ScalarTensor::Tags::ScalarMass>;
  using return_type = Scalar<DataType>;
  //   static constexpr void (*function)(const gsl::not_null<return_type*>
  //   result,
  //                                     const Scalar<DataType>&) =
  //       &compute_target_psi;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataType>&,
      const Scalar<DataType>&, const Scalar<DataType>&, const double,
      const double, const double) =
      &ScalarTensor::compute_scalar_curvature_source;
  using base = ScalarDriverSource;
};

}  // namespace fe::ScalarDriver::Tags
