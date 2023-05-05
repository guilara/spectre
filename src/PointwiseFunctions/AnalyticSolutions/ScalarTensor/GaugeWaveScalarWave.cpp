// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/GaugeWaveScalarWave.hpp"

#include <cmath>
#include <cstddef>

#include "DataStructures/DataVector.hpp"     // IWYU pragma: keep
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"

// IWYU pragma:  no_include "DataStructures/Tensor/TypeAliases.hpp"

namespace ScalarTensor::Solutions {

GaugeWaveScalarWave::GaugeWaveScalarWave(
    double amplitude_gauge_wave, double wavelength_gauge_wave,
    double amplitude, std::array<double, 3_st> wave_vector,
    std::array<double, 3_st> center,
    std::unique_ptr<MathFunction<1, Frame::Inertial>> profile)
    : amplitude_gauge_wave_(std::move(amplitude_gauge_wave)),
      wavelength_gauge_wave_(std::move(wavelength_gauge_wave)),
      amplitude_(std::move(amplitude)),
      wave_vector_(std::move(wave_vector)),
      center_(std::move(center)),
      omega_(magnitude(wave_vector_)),
      profile_(std::move(profile)),
      background_spacetime_(gr::Solutions::GaugeWave<3_st>{
                                                    amplitude_gauge_wave_,
                                                    wavelength_gauge_wave_}) {}

GaugeWaveScalarWave::GaugeWaveScalarWave(const GaugeWaveScalarWave& other)
    : evolution::initial_data::InitialData(other),
      amplitude_gauge_wave_(other.amplitude_gauge_wave_),
      wavelength_gauge_wave_(other.wavelength_gauge_wave_),
      amplitude_(other.amplitude_),
      wave_vector_(other.wave_vector_),
      center_(other.center_),
      omega_(magnitude(wave_vector_)),
      profile_(other.profile_->get_clone()),
      background_spacetime_(other.background_spacetime_) {}

GaugeWaveScalarWave& GaugeWaveScalarWave::operator=(
    const GaugeWaveScalarWave& other) {
  amplitude_gauge_wave_ = other.amplitude_gauge_wave_;
  wavelength_gauge_wave_ = other.wavelength_gauge_wave_;
  amplitude_ = other.amplitude_;
  wave_vector_ = other.wave_vector_;
  center_ = other.center_;
  omega_ = magnitude(wave_vector_);
  profile_ = other.profile_->get_clone();
  background_spacetime_ = background_spacetime_;
  return *this;
}

std::unique_ptr<evolution::initial_data::InitialData>
GaugeWaveScalarWave::get_clone() const {
  return std::make_unique<GaugeWaveScalarWave>(*this);
}

GaugeWaveScalarWave::GaugeWaveScalarWave(CkMigrateMessage* msg)
    : InitialData(msg) {}

void GaugeWaveScalarWave::pup(PUP::er& p) {
  InitialData::pup(p);
  p | amplitude_gauge_wave_;
  p | wavelength_gauge_wave_;
  p | amplitude_;
  p | wave_vector_;
  p | center_;
  p | profile_;
  p | omega_;
  p | background_spacetime_;
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Psi> GaugeWaveScalarWave::variables(
    const tnsr::I<DataType, 3>& x, double t,
    tmpl::list<CurvedScalarWave::Tags::Psi> /*meta*/) const {
    return {Scalar<DataType>(profile_->operator()(u(x, t)))};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Phi<3_st>>
GaugeWaveScalarWave::variables(
    const tnsr::I<DataType, 3>& x, double t,
    tmpl::list<CurvedScalarWave::Tags::Phi<3_st>> /*meta*/) const {
  auto result = make_with_value<tnsr::i<DataType, 3_st>>(x, 0.0);
  const auto du = profile_->first_deriv(u(x, t));
  for (size_t i = 0; i < 3_st; ++i) {
    result.get(i) = gsl::at(wave_vector_, i) * du;
    result.get(i) *= 2.0 * M_PI;
  }
  return {result};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Pi> GaugeWaveScalarWave::variables(
    const tnsr::I<DataType, 3>& x, double t,
    tmpl::list<CurvedScalarWave::Tags::Pi> /*meta*/) const {
  auto result = Scalar<DataType>(-omega_ * profile_->first_deriv(u(x, t)));
  result.get() *= 2.0 * M_PI;
  // Retrieve lapse from the backround spacetime
  auto spacetime_vars = background_spacetime_.variables(x, t,
                                  tmpl::list<gr::Tags::Lapse<DataType>>{});
  auto lapse = get<gr::Tags::Lapse<DataType>>(spacetime_vars);
  // For the Gauge Wave spacetime the shift is zero
  result.get() *= -1.0/lapse.get();
  return {result};
}

template <typename DataType>
DataType GaugeWaveScalarWave::u(const tnsr::I<DataType, 3_st>& x,
                                const double t) const {
  auto result = make_with_value<DataType>(x, -omega_ * t);
  for (size_t d = 0; d < 3_st; ++d) {
    result += gsl::at(wave_vector_, d) * (x.get(d) - gsl::at(center_, d));
  }
  result *= 2.0 * M_PI;
  return result;
}

PUP::able::PUP_ID GaugeWaveScalarWave::my_PUP_ID = 0;

bool operator==(const GaugeWaveScalarWave& lhs,
                const GaugeWaveScalarWave& rhs) {
  return (lhs.amplitude_gauge_wave_ == rhs.amplitude_gauge_wave_) and
         (lhs.wavelength_gauge_wave_ == rhs.wavelength_gauge_wave_) and
         (lhs.amplitude_ == rhs.amplitude_) and
         (lhs.wave_vector_ == rhs.wave_vector_) and
         (lhs.center_ == rhs.center_) and
         (*(lhs.profile_) == *(rhs.profile_)) and (lhs.omega_ == rhs.omega_);
}

bool operator!=(const GaugeWaveScalarWave& lhs,
                const GaugeWaveScalarWave& rhs) {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE_SCALARS(_, data)                                      \
  template tuples::TaggedTuple<TAG(data)> GaugeWaveScalarWave::variables( \
      const tnsr::I<DTYPE(data), 3>& x, double t,                         \
      tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_SCALARS, (DataVector),
                        (CurvedScalarWave::Tags::Psi,
                         CurvedScalarWave::Tags::Pi))

#define INSTANTIATE_VECTORS(_, data)                                      \
  template tuples::TaggedTuple<TAG(data)> GaugeWaveScalarWave::variables( \
      const tnsr::I<DTYPE(data), 3>& x, double t,                         \
      tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_VECTORS, (DataVector),
                        (CurvedScalarWave::Tags::Phi<3_st>))

template DataVector GaugeWaveScalarWave::u(                               \
      const tnsr::I<DataVector, 3_st>& x, const double t) const;

#undef DTYPE
#undef TAG
#undef INSTANTIATE_SCALARS
#undef INSTANTIATE_VECTORS
#undef INSTANTIATE
}  // namespace ScalarTensor::Solutions
