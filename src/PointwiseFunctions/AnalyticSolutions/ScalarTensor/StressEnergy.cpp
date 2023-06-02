// Distributed under the MIT License.
// See LICENSE.txt for details.

namespace ScalarTensor {
template <typename DataType>
void energy_density(gsl::not_null<Scalar<DataType>*> result,
                    const Scalar<DataType>& scalar_psi,
                    const tnsr::i<DataType, 3>& scalar_phi,
                    const Scalar<DataType>& scalar_pi) {
  //
}

template <typename DataType>
void momentum_density(gsl::not_null<tnsr::I<DataType, 3>*> result,
                      const Scalar<DataType>& scalar_psi,
                      const tnsr::i<DataType, 3>& scalar_phi,
                      const Scalar<DataType>& scalar_pi) {
  //
}

template <typename DataType>
void stress_trace(gsl::not_null<Scalar<DataType>*> result,
                  const Scalar<DataType>& scalar_psi,
                  const tnsr::i<DataType, 3>& scalar_phi,
                  const Scalar<DataType>& scalar_pi) {
  //
}
}  // namespace ScalarTensor
