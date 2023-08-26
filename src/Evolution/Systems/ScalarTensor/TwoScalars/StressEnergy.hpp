// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/TwoScalars/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor::TwoScalars {

void trace_reversed_stress_energy(
    gsl::not_null<tnsr::aa<DataVector, 3_st>*> stress_energy,
    const Scalar<DataVector>& pi_scalar,
    const tnsr::i<DataVector, 3_st>& phi_scalar,
    const Scalar<DataVector>& pi_scalar_2,
    const tnsr::i<DataVector, 3_st>& phi_scalar_2,
    const Scalar<DataVector>& lapse);

namespace Tags {

/*!
 * \brief Compute tag for the trace reversed stress energy tensor.
 *
 * \details Compute using ScalarTensor::trace_reversed_stress_energy.
 */
struct TraceReversedStressEnergyCompute
    : ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, 3_st,
                                                    Frame::Inertial>,
      db::ComputeTag {
  static constexpr size_t Dim = 3;
  using argument_tags = tmpl::list<Csw<CurvedScalarWave::Tags::Pi, 1>,
                                   Csw<CurvedScalarWave::Tags::Phi<Dim>, 1>,
                                   Csw<CurvedScalarWave::Tags::Pi, 2>,
                                   Csw<CurvedScalarWave::Tags::Phi<Dim>, 2>,
                                   gr::Tags::Lapse<DataVector>>;
  using return_type = tnsr::aa<DataVector, Dim, Frame::Inertial>;
  static constexpr void (*function)(
      const gsl::not_null<tnsr::aa<DataVector, Dim>*> result,
      const Scalar<DataVector>&, const tnsr::i<DataVector, Dim>&,
      const Scalar<DataVector>&, const tnsr::i<DataVector, Dim>&,
      const Scalar<DataVector>&) = &trace_reversed_stress_energy;
  using base = ScalarTensor::Tags::TraceReversedStressEnergy<DataVector, Dim,
                                                             Frame::Inertial>;
};
}  // namespace Tags

}  // namespace ScalarTensor::TwoScalars
