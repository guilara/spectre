// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

// Here we write the functions computing the scalar stress-energy tensor
// projections for the initial data solves (i.e. energy density, momentum
// density, and trace of the spatial stress-tensor

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarTensor {
template <typename DataType>
void energy_density(gsl::not_null<Scalar<DataType>*> result,
                    const Scalar<DataType>& scalar_psi,
                    const tnsr::i<DataType, 3>& scalar_phi,
                    const Scalar<DataType>& scalar_pi);

template <typename DataType>
void momentum_density(gsl::not_null<tnsr::I<DataType, 3>*> result,
                      const Scalar<DataType>& scalar_psi,
                      const tnsr::i<DataType, 3>& scalar_phi,
                      const Scalar<DataType>& scalar_pi);

template <typename DataType>
void stress_trace(gsl::not_null<Scalar<DataType>*> result,
                  const Scalar<DataType>& scalar_psi,
                  const tnsr::i<DataType, 3>& scalar_phi,
                  const Scalar<DataType>& scalar_pi);
}  // namespace ScalarTensor
