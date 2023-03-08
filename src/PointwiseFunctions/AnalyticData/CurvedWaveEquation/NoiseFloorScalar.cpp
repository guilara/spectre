// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/CurvedWaveEquation/NoiseFloorScalar.hpp"

#include <complex>
#include <cstddef>
#include <pup.h>
#include <random>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
//
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
//
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace CurvedScalarWave::AnalyticData {

NoiseFloorScalar::NoiseFloorScalar(const double amplitude,
                                             const Options::Context& context)
    : amplitude_(amplitude) {
  if (amplitude_ <= 0.) {
    PARSE_ERROR(context,
                "The amplitude must be greater than 0 but is " << amplitude_);
  }
}

tuples::TaggedTuple<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
                    CurvedScalarWave::Tags::Phi<3>>
NoiseFloorScalar::variables(const tnsr::I<DataVector, 3>& x,
                                 tags /*meta*/) const {
  // Need to figure out this
  // std::mt19937 generator{111};
  // std::uniform_real_distribution<> distribution(-0.5, 0.5);
  // Scalar<DataVector> pi = make_with_random_values<Scalar<DataVector>>(
  //     make_not_null(&generator), make_not_null(&distribution), x);

  // Make with constant value
  Scalar<DataVector> psi{amplitude_};
  Scalar<DataVector> pi{amplitude_};
  tnsr::i<DataVector, 3> phi{amplitude_};

  return tuples::TaggedTuple<CurvedScalarWave::Tags::Psi,
                             CurvedScalarWave::Tags::Pi,
                             CurvedScalarWave::Tags::Phi<3>>{
      std::move(psi), std::move(pi), std::move(phi)};
}

void NoiseFloorScalar::pup(PUP::er& p) {
  p | amplitude_;
}

bool operator==(const NoiseFloorScalar& lhs,
                const NoiseFloorScalar& rhs) {
  return lhs.amplitude_ == rhs.amplitude_;
}
bool operator!=(const NoiseFloorScalar& lhs,
                const NoiseFloorScalar& rhs) {
  return not(lhs == rhs);
}
}  // namespace CurvedScalarWave::AnalyticData
