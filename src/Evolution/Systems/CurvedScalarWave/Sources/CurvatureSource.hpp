// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Sources/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "Utilities/Gsl.hpp"

namespace CurvedScalarWave::Sources {

/*!
 * \brief Compute the coupling function entering the scalar source term for the
 * CurvedScalarWave system.
 *
 * \details The scalar source term depends on the problem at hand.
 * Here we write a scalar source given by the curvature of the background
 * spacetime and depending on two coupling parameters.
 *
 * \f[
 * \mathrm{scalar source} = f'(\psi) E_{ab} E^{ab},
 * \f]
 *
 * where \f$ f'(\psi) = p_1 \psi + p_2 \psi^2 \f$.
 */
Scalar<DataVector> coupling_function_prime(const Scalar<DataVector>& psi,
                                            const double first_coupling_psi,
                                            const double second_coupling_psi);

/*!
 * \brief Compute the scalar source term for the CurvedScalarWave system.
 *
 * \details The scalar source term depends on the problem at hand.
 * Here we write a scalar source given by the curvature of the background
 * spacetime,
 *
 * \f[
 * \mathrm{scalar source} = f'(\psi) E_{ab} E^{ab},
 * \f]
 *
 * where \f$ f'(\psi) = p_1 \psi + p_2 \psi^2 \f$.
 */
void compute_scalar_curvature_source(
    gsl::not_null<Scalar<DataVector>*> scalar_source,
    const Scalar<DataVector>& weyl_electric_scalar,
    const Scalar<DataVector>& psi, const double first_coupling_psi,
    const double second_coupling_psi);

namespace Tags {

/*!
 * \brief Compute tag for the scalar source given by the background curvature.
 *
 * \details Call compute_scalar_curvature_source. Needs that WeylElectric is in
 * data box.
 */
template <size_t SpatialDim, typename Frame, typename DataType>
struct ScalarCurvatureSourceCompute : ScalarSource, db::ComputeTag {
//   using argument_tags = tmpl::list<
//       gr::Tags::WeylElectricScalarCompute<SpatialDim, Frame, DataType>,
//       CurvedScalarWave::Sources::Tags::ScalarMass>;
using argument_tags =
    tmpl::list<gr::Tags::WeylElectricScalarCompute<SpatialDim, Frame, DataType>,
               CurvedScalarWave::Tags::Psi,
               CurvedScalarWave::Sources::Tags::ScalarFirstCouplingParameter,
               CurvedScalarWave::Sources::Tags::ScalarSecondCouplingParameter>;
using return_type = Scalar<DataVector>;
static constexpr void (*function)(
    const gsl::not_null<return_type*> result, const Scalar<DataVector>&,
    const Scalar<DataVector>&, const double,
    const double) = &compute_scalar_curvature_source;
using base = ScalarSource;
};

}  // namespace Tags

}  // namespace CurvedScalarWave::Sources
