// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/BoundaryCorrections/UpwindPenalty.hpp"

#include <memory>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Characteristics.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"  // IWYU pragma: keep

namespace fe::ScalarDriver::BoundaryCorrections {
UpwindPenalty::UpwindPenalty(CkMigrateMessage* msg) : BoundaryCorrection(msg) {}

std::unique_ptr<BoundaryCorrection> UpwindPenalty::get_clone() const {
  return std::make_unique<UpwindPenalty>(*this);
}

void UpwindPenalty::pup(PUP::er& p) {
  BoundaryCorrection::pup(p);
  p | boundary_correction_for_scalar_;
}

// NOLINTNEXTLINE
PUP::able::PUP_ID UpwindPenalty::my_PUP_ID = 0;

}  // namespace fe::ScalarDriver::BoundaryCorrections
