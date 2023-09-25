// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/ScalarTensor/SphericalKerrSchildSH.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include "DataStructures/DataVector.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/Spectral/SwshInterpolation.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

// IWYU pragma:  no_include "DataStructures/Tensor/TypeAliases.hpp"

namespace ScalarTensor::AnalyticData {

SphericalKerrSchildSH::SphericalKerrSchildSH(
    const double mass, const std::array<double, 3>& dimensionless_spin,
    const double amplitude, const double radius, const double width,
    const std::pair<size_t, int> mode) {
  mass_ = mass;
  dimensionless_spin_ = dimensionless_spin;
  amplitude_ = amplitude;
  radius_ = radius;
  width_sq_ = square(width);
  mode_ = mode;
  background_spacetime_ = gr::Solutions::SphericalKerrSchild{
      mass_, dimensionless_spin_,
      // Center
      std::array<double, 3>{{0.0, 0.0, 0.0}}};
}

std::unique_ptr<evolution::initial_data::InitialData>
SphericalKerrSchildSH::get_clone() const {
  return std::make_unique<SphericalKerrSchildSH>(*this);
}

SphericalKerrSchildSH::SphericalKerrSchildSH(CkMigrateMessage* msg)
    : InitialData(msg) {}

void SphericalKerrSchildSH::pup(PUP::er& p) {
  InitialData::pup(p);
  p | mass_;
  p | dimensionless_spin_;
  p | amplitude_;
  p | radius_;
  p | width_sq_;
  p | mode_;
  p | background_spacetime_;
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Psi>
SphericalKerrSchildSH::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Psi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Phi<3_st>>
SphericalKerrSchildSH::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Phi<3_st>> /*meta*/) const {
  return {make_with_value<tnsr::i<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Pi>
SphericalKerrSchildSH::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Pi> /*meta*/) const {
  Scalar<DataType> pi = make_with_value<Scalar<DataType>>(x, 0.0);
  get(pi) += get(magnitude(x)) - radius_;
  get(pi) = exp(-get(pi) * get(pi) / width_sq_);
  const Spectral::Swsh::SpinWeightedSphericalHarmonic spherical_harmonic(
      0, mode_.first, mode_.second);
  const auto theta = atan2(hypot(x[0], x[1]), x[2]);
  const auto phi = atan2(x[1], x[0]);
  get(pi) *= real(spherical_harmonic.evaluate(theta, phi, sin(theta / 2.),
                                              cos(theta / 2.)));
  get(pi) *= amplitude_;
  return pi;
}

PUP::able::PUP_ID SphericalKerrSchildSH::my_PUP_ID = 0;

bool operator==(const SphericalKerrSchildSH& lhs,
                const SphericalKerrSchildSH& rhs) {
  return lhs.amplitude_ == rhs.amplitude_ and
         lhs.dimensionless_spin_ == rhs.dimensionless_spin_ and
         lhs.radius_ == rhs.radius_ and lhs.width_sq_ == rhs.width_sq_ and
         lhs.mode_ == rhs.mode_ and lhs.mass_ == rhs.mass_ and
         lhs.background_spacetime_ == rhs.background_spacetime_;
}

bool operator!=(const SphericalKerrSchildSH& lhs,
                const SphericalKerrSchildSH& rhs) {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE_SCALARS(_, data)                                        \
  template tuples::TaggedTuple<TAG(data)> SphericalKerrSchildSH::variables( \
      const tnsr::I<DTYPE(data), 3>& x, tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_SCALARS, (DataVector),
                        (CurvedScalarWave::Tags::Psi,
                         CurvedScalarWave::Tags::Pi))

#define INSTANTIATE_VECTORS(_, data)                                        \
  template tuples::TaggedTuple<TAG(data)> SphericalKerrSchildSH::variables( \
      const tnsr::I<DTYPE(data), 3>& x, tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_VECTORS, (DataVector),
                        (CurvedScalarWave::Tags::Phi<3_st>))

#undef DTYPE
#undef TAG
#undef INSTANTIATE_SCALARS
#undef INSTANTIATE_VECTORS
}  // namespace ScalarTensor::AnalyticData
