// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor::Sources {

void add_scalar_source_to_dt_pi_scalar(gsl::not_null<Scalar<DataVector>*> dt_pi,
                                       const Scalar<DataVector>& scalar_source,
                                       const Scalar<DataVector>& lapse);

void compute_scalar_source(gsl::not_null<Scalar<DataVector>*> scalar_source,
                           const Scalar<DataVector>& psi,
                           const double mass_psi);

void compute_scalar_curvature_source(
    gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi);

Scalar<DataVector> compute_scalar_curvature_source(
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& weyl_magnetic_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi);

void multiply_by_coupling_function_prime_quartic(
    gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi);

namespace Tags {

/*!
 * \brief Compute tag for the scalar source.
 *
 * \details Call compute_scalar_source.
 */
struct ScalarSourceCompute : ScalarSource, db::ComputeTag {
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Psi,
                                   ScalarTensor::Sources::Tags::ScalarMass>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(const gsl::not_null<return_type*> result,
                                    const Scalar<DataVector>&,
                                    const double) = &compute_scalar_source;
  using base = ScalarSource;
};

/*!
 * \brief Compute tag for the scalar source given by the background curvature.
 *
 * \details Call compute_scalar_curvature_source. Needs that WeylElectric is in
 * data box.
 */
template <size_t SpatialDim, typename Frame, typename DataType>
struct ScalarCurvatureSourceCompute : ScalarSource, db::ComputeTag {
  using argument_tags = tmpl::list<
      gr::Tags::WeylElectricScalar<DataType>,
      gr::Tags::WeylMagneticScalar<DataType>,
      CurvedScalarWave::Tags::Psi,
      ScalarTensor::Sources::Tags::ScalarFirstCouplingParameter,
      ScalarTensor::Sources::Tags::ScalarSecondCouplingParameter,
      ScalarTensor::Sources::Tags::ScalarMass>;
  using return_type = Scalar<DataVector>;
  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataVector>&,
      const Scalar<DataVector>&, const Scalar<DataVector>&, const double,
      const double, const double) = &compute_scalar_curvature_source;
  using base = ScalarSource;
};

}  // namespace Tags

}  // namespace ScalarTensor::Sources
