// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/MinkowskiScalarWave.hpp"

#include <cmath>
#include <cstddef>

#include "DataStructures/DataVector.hpp"     // IWYU pragma: keep
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"

// IWYU pragma:  no_include "DataStructures/Tensor/TypeAliases.hpp"

namespace ScalarTensor::Solutions {

MinkowskiScalarWave::MinkowskiScalarWave(
    double amplitude, std::array<double, 3_st> wave_vector,
    std::array<double, 3_st> center,
    std::unique_ptr<MathFunction<1, Frame::Inertial>> profile)
    : amplitude_(std::move(amplitude)),
      wave_vector_(std::move(wave_vector)),
      center_(std::move(center)),
      omega_(magnitude(wave_vector_)),
      profile_(std::move(profile)) {}

MinkowskiScalarWave::MinkowskiScalarWave(const MinkowskiScalarWave& other)
    : evolution::initial_data::InitialData(other),
      amplitude_(other.amplitude_),
      wave_vector_(other.wave_vector_),
      center_(other.center_),
      omega_(magnitude(wave_vector_)),
      profile_(other.profile_->get_clone()) {}

MinkowskiScalarWave& MinkowskiScalarWave::operator=(
    const MinkowskiScalarWave& other) {
  amplitude_ = other.amplitude_;
  wave_vector_ = other.wave_vector_;
  center_ = other.center_;
  omega_ = magnitude(wave_vector_);
  profile_ = other.profile_->get_clone();
  return *this;
}

std::unique_ptr<evolution::initial_data::InitialData>
MinkowskiScalarWave::get_clone() const {
  return std::make_unique<MinkowskiScalarWave>(*this);
}

MinkowskiScalarWave::MinkowskiScalarWave(CkMigrateMessage* msg)
    : InitialData(msg) {}

void MinkowskiScalarWave::pup(PUP::er& p) {
  InitialData::pup(p);
  p | amplitude_;
  p | wave_vector_;
  p | center_;
  p | profile_;
  p | background_spacetime_;
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Psi> MinkowskiScalarWave::variables(
    const tnsr::I<DataType, 3>& x, double t,
    tmpl::list<CurvedScalarWave::Tags::Psi> /*meta*/) const {
    return {Scalar<DataType>(profile_->operator()(u(x, t)))};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Phi<3_st>>
MinkowskiScalarWave::variables(
    const tnsr::I<DataType, 3>& x, double t,
    tmpl::list<CurvedScalarWave::Tags::Phi<3_st>> /*meta*/) const {
  auto result = make_with_value<tnsr::i<DataType, 3_st>>(x, 0.0);
  const auto du = profile_->first_deriv(u(x, t));
  for (size_t i = 0; i < 3_st; ++i) {
    result.get(i) = gsl::at(wave_vector_, i) * du;
  }
  return {result};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Pi> MinkowskiScalarWave::variables(
    const tnsr::I<DataType, 3>& x, double t,
    tmpl::list<CurvedScalarWave::Tags::Pi> /*meta*/) const {
  auto result = Scalar<DataType>(-omega_ * profile_->first_deriv(u(x, t)));
  result.get() *= -1.0;
  return {result};
}

template <typename DataType>
DataType MinkowskiScalarWave::u(const tnsr::I<DataType, 3_st>& x,
                                const double t) const {
  auto result = make_with_value<DataType>(x, -omega_ * t);
  for (size_t d = 0; d < 3_st; ++d) {
    result += gsl::at(wave_vector_, d) * (x.get(d) - gsl::at(center_, d));
  }
  return result;
}

PUP::able::PUP_ID MinkowskiScalarWave::my_PUP_ID = 0;

bool operator==(const MinkowskiScalarWave& lhs,
                const MinkowskiScalarWave& rhs) {
  return (lhs.amplitude_ == rhs.amplitude_) and
         (lhs.wave_vector_ == rhs.wave_vector_) and
         (lhs.center_ == rhs.center_) and
         (*(lhs.profile_) == *(rhs.profile_)) and (lhs.omega_ == rhs.omega_);
}

bool operator!=(const MinkowskiScalarWave& lhs,
                const MinkowskiScalarWave& rhs) {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE_SCALARS(_, data)                                      \
  template tuples::TaggedTuple<TAG(data)> MinkowskiScalarWave::variables( \
      const tnsr::I<DTYPE(data), 3>& x, double t,                         \
      tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_SCALARS, (DataVector),
                        (CurvedScalarWave::Tags::Psi,
                         CurvedScalarWave::Tags::Pi))

#define INSTANTIATE_VECTORS(_, data)                                      \
  template tuples::TaggedTuple<TAG(data)> MinkowskiScalarWave::variables( \
      const tnsr::I<DTYPE(data), 3>& x, double t,                         \
      tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_VECTORS, (DataVector),
                        (CurvedScalarWave::Tags::Phi<3_st>))

template DataVector MinkowskiScalarWave::u(                               \
      const tnsr::I<DataVector, 3_st>& x, const double t) const;

#undef DTYPE
#undef TAG
#undef INSTANTIATE_SCALARS
#undef INSTANTIATE_VECTORS
#undef INSTANTIATE
}  // namespace ScalarTensor::Solutions
