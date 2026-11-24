// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <random>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"
#include "PointwiseFunctions/SpecialRelativity/LorentzBoostMatrix.hpp"
#include "PointwiseFunctions/SpecialRelativity/LorentzContractionMatrix.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/StdHelpers.hpp"

namespace {
template <size_t SpatialDim>
tnsr::Ij<double, SpatialDim, Frame::NoFrame> identity_matrix() {
  tnsr::Ij<double, SpatialDim, Frame::NoFrame> identity{0.};
  for (size_t i = 0; i < SpatialDim; ++i) {
    identity.get(i, i) = 1.;
  }
  return identity;
}

template <size_t SpatialDim>
std::array<double, SpatialDim> matrix_times_vector(
    const tnsr::Ij<double, SpatialDim, Frame::NoFrame>& matrix,
    const std::array<double, SpatialDim>& vector) {
  std::array<double, SpatialDim> result = make_array<SpatialDim>(0.);
  for (size_t i = 0; i < SpatialDim; ++i) {
    for (size_t j = 0; j < SpatialDim; ++j) {
      gsl::at(result, i) += matrix.get(i, j) * gsl::at(vector, j);
    }
  }
  return result;
}

template <size_t SpatialDim>
void test_zero_velocity() {
  const auto velocity = make_array<SpatialDim>(0.);
  CHECK(sr::lorentz_contraction_matrix(velocity) ==
        identity_matrix<SpatialDim>());
  CHECK(sr::inverse_lorentz_contraction_matrix(velocity) ==
        identity_matrix<SpatialDim>());
}

template <size_t SpatialDim>
void test_random_velocity(const gsl::not_null<std::mt19937*> generator) {
  std::uniform_real_distribution<> distribution(-1., 1.);
  std::uniform_real_distribution<> speed_distribution(0., 0.99);
  auto velocity = make_with_random_values<std::array<double, SpatialDim>>(
      generator, make_not_null(&distribution));
  double norm = 0.;
  for (size_t i = 0; i < SpatialDim; ++i) {
    norm += square(gsl::at(velocity, i));
  }
  norm = sqrt(norm);
  const double speed = speed_distribution(*generator);
  for (size_t i = 0; i < SpatialDim; ++i) {
    gsl::at(velocity, i) *= speed / norm;
  }
  const double lorentz_factor = 1. / sqrt(1. - square(speed));
  CAPTURE(velocity);

  const auto contraction_matrix = sr::lorentz_contraction_matrix(velocity);
  const auto inverse_contraction_matrix =
      sr::inverse_lorentz_contraction_matrix(velocity);

  // Matches the spatial block of the Lorentz boost matrix
  const auto boost_matrix = sr::lorentz_boost_matrix(velocity);
  for (size_t i = 0; i < SpatialDim; ++i) {
    for (size_t j = 0; j < SpatialDim; ++j) {
      CHECK(contraction_matrix.get(i, j) ==
            approx(boost_matrix.get(i + 1, j + 1)));
    }
  }

  // Product with the inverse is the identity
  tnsr::Ij<double, SpatialDim, Frame::NoFrame> product{0.};
  for (size_t i = 0; i < SpatialDim; ++i) {
    for (size_t j = 0; j < SpatialDim; ++j) {
      for (size_t k = 0; k < SpatialDim; ++k) {
        product.get(i, j) +=
            contraction_matrix.get(i, k) * inverse_contraction_matrix.get(k, j);
      }
    }
  }
  CHECK_ITERABLE_APPROX(product, identity_matrix<SpatialDim>());

  // Stretches (contracts) along the velocity by gamma (1/gamma)
  auto expected = velocity;
  for (size_t i = 0; i < SpatialDim; ++i) {
    gsl::at(expected, i) *= lorentz_factor;
  }
  CHECK_ITERABLE_APPROX(matrix_times_vector(contraction_matrix, velocity),
                        expected);
  for (size_t i = 0; i < SpatialDim; ++i) {
    gsl::at(expected, i) = gsl::at(velocity, i) / lorentz_factor;
  }
  CHECK_ITERABLE_APPROX(
      matrix_times_vector(inverse_contraction_matrix, velocity), expected);

  // Leaves directions orthogonal to the velocity unchanged
  if constexpr (SpatialDim > 1) {
    auto orthogonal = make_array<SpatialDim>(0.);
    orthogonal[0] = velocity[1];
    orthogonal[1] = -velocity[0];
    CHECK_ITERABLE_APPROX(matrix_times_vector(contraction_matrix, orthogonal),
                          orthogonal);
    CHECK_ITERABLE_APPROX(
        matrix_times_vector(inverse_contraction_matrix, orthogonal),
        orthogonal);
  }
}

void test_small_velocity() {
  // For v << 1 the off-diagonal elements are v_i v_j (1/2 + O(v^2)) for the
  // contraction matrix and -v_i v_j (1/2 + O(v^2)) for its inverse. They must
  // be accurate to full relative precision, which fails if they are computed
  // as differences of O(1) numbers.
  const double small = 1.e-9;
  const std::array<double, 3> velocity{{small, 2. * small, -3. * small}};
  const auto contraction_matrix = sr::lorentz_contraction_matrix(velocity);
  const auto inverse_contraction_matrix =
      sr::inverse_lorentz_contraction_matrix(velocity);
  Approx custom_approx = Approx::custom().epsilon(1.e-14).scale(0.);
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      if (i == j) {
        continue;
      }
      const double expected = 0.5 * gsl::at(velocity, i) * gsl::at(velocity, j);
      CHECK(contraction_matrix.get(i, j) == custom_approx(expected));
      CHECK(inverse_contraction_matrix.get(i, j) == custom_approx(-expected));
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.SpecialRelativity.LorentzContraction",
    "[Unit][PointwiseFunctions]") {
  MAKE_GENERATOR(generator);
  test_zero_velocity<1>();
  test_zero_velocity<2>();
  test_zero_velocity<3>();
  for (size_t i = 0; i < 100; ++i) {
    test_random_velocity<1>(make_not_null(&generator));
    test_random_velocity<2>(make_not_null(&generator));
    test_random_velocity<3>(make_not_null(&generator));
  }
  test_small_velocity();
}
