// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/SpecialRelativity/LorentzContractionMatrix.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/SpecialRelativity/LorentzBoostMatrix.hpp"
#include "Utilities/GenerateInstantiations.hpp"

namespace sr {
template <size_t SpatialDim>
tnsr::Ij<double, SpatialDim, Frame::NoFrame> lorentz_contraction_matrix(
    const std::array<double, SpatialDim>& velocity) {
  // The spatial block of the boost matrix is
  // delta_ij + gamma^2 / (1 + gamma) v_i v_j, computed without dividing by v^2
  const auto boost_matrix = lorentz_boost_matrix(velocity);
  tnsr::Ij<double, SpatialDim, Frame::NoFrame> contraction_matrix{};
  for (size_t i = 0; i < SpatialDim; ++i) {
    for (size_t j = 0; j < SpatialDim; ++j) {
      contraction_matrix.get(i, j) = boost_matrix.get(i + 1, j + 1);
    }
  }
  return contraction_matrix;
}

template <size_t SpatialDim>
tnsr::Ij<double, SpatialDim, Frame::NoFrame> inverse_lorentz_contraction_matrix(
    const std::array<double, SpatialDim>& velocity) {
  // (1 - 1/gamma) / v^2 = gamma / (1 + gamma). Writing the outer product with
  // the boost-matrix components gamma v_i avoids both the division by v^2 and
  // the cancellation in 1 - 1/gamma for small velocities. Subtracting the
  // identity from the spatial block of the boost matrix would instead lose
  // the significant digits of the O(v^2) term.
  const auto boost_matrix = lorentz_boost_matrix(velocity);
  const double lorentz_factor = get<0, 0>(boost_matrix);
  const double prefactor = 1. / (lorentz_factor * (1. + lorentz_factor));
  tnsr::Ij<double, SpatialDim, Frame::NoFrame> inverse_contraction_matrix{};
  for (size_t i = 0; i < SpatialDim; ++i) {
    for (size_t j = 0; j < SpatialDim; ++j) {
      inverse_contraction_matrix.get(i, j) =
          -prefactor * boost_matrix.get(i + 1, 0) * boost_matrix.get(j + 1, 0);
    }
    inverse_contraction_matrix.get(i, i) += 1.;
  }
  return inverse_contraction_matrix;
}
}  // namespace sr

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                               \
  template tnsr::Ij<double, DIM(data), Frame::NoFrame>                     \
  sr::lorentz_contraction_matrix(                                          \
      const std::array<double, DIM(data)>& velocity);                      \
  template tnsr::Ij<double, DIM(data), Frame::NoFrame>                     \
  sr::inverse_lorentz_contraction_matrix(                                  \
      const std::array<double, DIM(data)>& velocity);

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef DIM
#undef INSTANTIATE
