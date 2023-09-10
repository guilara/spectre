// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Actions/InitializeConstraintGammas.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::ScalarDriver::Initialization {

void InitializeConstraintDampingGammasGaussian::apply(
    const gsl::not_null<Scalar<DataVector>*> gamma1,
    const gsl::not_null<Scalar<DataVector>*> gamma2,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const double amplitude_gaussian, const double sigma_gaussian,
    const double offset_gaussian) {
  const size_t number_of_grid_points = get<0>(inertial_coords).size();
  *gamma1 = Scalar<DataVector>{number_of_grid_points, 0.};
  const double beta_gaussian = 1.0 / sigma_gaussian;
  const auto radius = magnitude(inertial_coords);
  get(*gamma2) =
      amplitude_gaussian * exp(-square(beta_gaussian * radius.get())) +
      offset_gaussian;
}

}  // namespace fe::ScalarDriver::Initialization
