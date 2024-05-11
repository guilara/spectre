// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/FixedScalarTensor/FixedSGB/KerrSphericalHarmonic.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <cstddef>

#include "DataStructures/DataVector.hpp"  // IWYU pragma: keep
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "NumericalAlgorithms/Spectral/SwshInterpolation.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

// IWYU pragma:  no_include "DataStructures/Tensor/TypeAliases.hpp"

namespace fe::DecoupledScalar::AnalyticData {

KerrSphericalHarmonic::KerrSphericalHarmonic(
    const double mass, const std::array<double, 3>& dimensionless_spin,
    const std::array<double, 3>& center,
    const std::array<double, 3>& boost_velocity, const double amplitude,
    const double amplitude_of_driver, const double radius, const double width,
    const std::pair<size_t, int> mode) {
  mass_ = mass;
  dimensionless_spin_ = dimensionless_spin;
  center_ = center;
  boost_velocity_ = boost_velocity;
  amplitude_ = amplitude;
  amplitude_of_driver_ = amplitude_of_driver;
  radius_ = radius;
  width_sq_ = square(width);
  mode_ = mode;
  background_spacetime_ = gr::Solutions::KerrSchild{mass_, dimensionless_spin_,
                                                    center_, boost_velocity_};
}

std::unique_ptr<evolution::initial_data::InitialData>
KerrSphericalHarmonic::get_clone() const {
  return std::make_unique<KerrSphericalHarmonic>(*this);
}

KerrSphericalHarmonic::KerrSphericalHarmonic(CkMigrateMessage* msg)
    : InitialData(msg) {}

void KerrSphericalHarmonic::pup(PUP::er& p) {
  InitialData::pup(p);
  p | mass_;
  p | dimensionless_spin_;
  p | center_;
  p | boost_velocity_;
  p | amplitude_;
  p | amplitude_of_driver_;
  p | radius_;
  p | width_sq_;
  p | mode_;
  p | background_spacetime_;
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Psi>
KerrSphericalHarmonic::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Psi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Phi<3_st>>
KerrSphericalHarmonic::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Phi<3_st>> /*meta*/) const {
  return {make_with_value<tnsr::i<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Pi>
KerrSphericalHarmonic::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Pi> /*meta*/) const {
  Scalar<DataType> pi = make_with_value<Scalar<DataType>>(x, 0.0);
  get(pi) += get(magnitude(x)) - radius_;
  get(pi) = exp(-square(get(pi)) / width_sq_);
  const Spectral::Swsh::SpinWeightedSphericalHarmonic spherical_harmonic(
      0, mode_.first, mode_.second);
  const auto theta = atan2(hypot(x[0], x[1]), x[2]);
  const auto phi = atan2(x[1], x[0]);
  get(pi) *= real(spherical_harmonic.evaluate(theta, phi, sin(0.5 * theta),
                                              cos(0.5 * theta)));
  get(pi) *= amplitude_;
  return pi;
}

// Scalar driver
template <typename DataType>
tuples::TaggedTuple<fe::ScalarTensorDriver::Tags::Psi>
KerrSphericalHarmonic::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<fe::ScalarTensorDriver::Tags::Psi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<fe::ScalarTensorDriver::Tags::TensorDriver<DataType, 3_st>>
KerrSphericalHarmonic::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<
        fe::ScalarTensorDriver::Tags::TensorDriver<DataType, 3_st>> /*meta*/)
    const {
  return {make_with_value<tnsr::aa<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<fe::ScalarTensorDriver::Tags::Pi<DataType, 3_st>>
KerrSphericalHarmonic::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<fe::ScalarTensorDriver::Tags::Pi<DataType, 3_st>> /*meta*/)
    const {
  return {make_with_value<tnsr::aa<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<fe::ScalarTensorDriver::Tags::PiScalar>
KerrSphericalHarmonic::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<fe::ScalarTensorDriver::Tags::PiScalar> /*meta*/) const {
  Scalar<DataType> pi = make_with_value<Scalar<DataType>>(x, 0.0);
  get(pi) += get(magnitude(x)) - radius_;
  get(pi) = exp(-square(get(pi)) / width_sq_);
  const Spectral::Swsh::SpinWeightedSphericalHarmonic spherical_harmonic(
      0, mode_.first, mode_.second);
  const auto theta = atan2(hypot(x[0], x[1]), x[2]);
  const auto phi = atan2(x[1], x[0]);
  get(pi) *= real(spherical_harmonic.evaluate(theta, phi, sin(0.5 * theta),
                                              cos(0.5 * theta)));
  get(pi) *= amplitude_of_driver_;
  return pi;
}

PUP::able::PUP_ID KerrSphericalHarmonic::my_PUP_ID = 0;

bool operator==(const KerrSphericalHarmonic& lhs,
                const KerrSphericalHarmonic& rhs) {
  return lhs.amplitude_ == rhs.amplitude_ and
         lhs.amplitude_of_driver_ == rhs.amplitude_of_driver_ and
         lhs.dimensionless_spin_ == rhs.dimensionless_spin_ and
         lhs.center_ == rhs.center_ and
         lhs.boost_velocity_ == rhs.boost_velocity_ and
         lhs.radius_ == rhs.radius_ and lhs.width_sq_ == rhs.width_sq_ and
         lhs.mode_ == rhs.mode_ and lhs.mass_ == rhs.mass_ and
         lhs.background_spacetime_ == rhs.background_spacetime_;
}

bool operator!=(const KerrSphericalHarmonic& lhs,
                const KerrSphericalHarmonic& rhs) {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE_SCALARS(_, data)                                        \
  template tuples::TaggedTuple<TAG(data)> KerrSphericalHarmonic::variables( \
      const tnsr::I<DTYPE(data), 3>& x, tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_SCALARS, (DataVector),
                        (CurvedScalarWave::Tags::Psi,
                         CurvedScalarWave::Tags::Pi,
                         fe::ScalarTensorDriver::Tags::Psi,
                         fe::ScalarTensorDriver::Tags::PiScalar))

#define INSTANTIATE_VECTORS(_, data)                                        \
  template tuples::TaggedTuple<TAG(data)> KerrSphericalHarmonic::variables( \
      const tnsr::I<DTYPE(data), 3>& x, tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_VECTORS, (DataVector),
                        (CurvedScalarWave::Tags::Phi<3_st>))

#define INSTANTIATE_TENSORS(_, data)                      \
  template tuples::TaggedTuple < TAG(data) < DTYPE(data), \
      3 >> KerrSphericalHarmonic::variables(              \
               const tnsr::I<DTYPE(data), 3>& x,          \
               tmpl::list < TAG(data) < DTYPE(data), 3 >> \
               /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_TENSORS, (DataVector),
                        (fe::ScalarTensorDriver::Tags::TensorDriver,
                         fe::ScalarTensorDriver::Tags::Pi))

#undef DTYPE
#undef TAG
#undef INSTANTIATE_SCALARS
#undef INSTANTIATE_VECTORS
#undef INSTANTIATE_TENSORS
}  // namespace fe::DecoupledScalar::AnalyticData
