// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/GaugeWaveConstantScalar.hpp"

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

// IWYU pragma:  no_include "DataStructures/Tensor/TypeAliases.hpp"

namespace ScalarTensor::Solutions {

GaugeWaveConstantScalar::GaugeWaveConstantScalar(double amplitude_gauge_wave,
                                                 double wavelength_gauge_wave,
                                                 double amplitude_scalar) {
  amplitude_gauge_wave_ = amplitude_gauge_wave;
  wavelength_gauge_wave_ = wavelength_gauge_wave;
  amplitude_scalar_ = amplitude_scalar;
  background_spacetime_ = gr::Solutions::GaugeWave<3_st>{amplitude_gauge_wave,
                                                         wavelength_gauge_wave};
}

std::unique_ptr<evolution::initial_data::InitialData>
GaugeWaveConstantScalar::get_clone() const {
  return std::make_unique<GaugeWaveConstantScalar>(*this);
}

GaugeWaveConstantScalar::GaugeWaveConstantScalar(CkMigrateMessage* msg)
    : InitialData(msg) {}

void GaugeWaveConstantScalar::pup(PUP::er& p) {
  InitialData::pup(p);
  p | amplitude_gauge_wave_;
  p | wavelength_gauge_wave_;
  p | amplitude_scalar_;
  p | background_spacetime_;
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Psi>
GaugeWaveConstantScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<CurvedScalarWave::Tags::Psi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, amplitude_scalar_)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Phi<3_st>>
GaugeWaveConstantScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<CurvedScalarWave::Tags::Phi<3_st>> /*meta*/) const {
  return {make_with_value<tnsr::i<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Pi>
GaugeWaveConstantScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<CurvedScalarWave::Tags::Pi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, 0.0)};
}

PUP::able::PUP_ID GaugeWaveConstantScalar::my_PUP_ID = 0;

bool operator==(const GaugeWaveConstantScalar& lhs,
                const GaugeWaveConstantScalar& rhs) {
  return lhs.background_spacetime_ == rhs.background_spacetime_;
}

bool operator!=(const GaugeWaveConstantScalar& lhs,
                const GaugeWaveConstantScalar& rhs) {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE_SCALARS(_, data)                                          \
  template tuples::TaggedTuple<TAG(data)> GaugeWaveConstantScalar::variables( \
      const tnsr::I<DTYPE(data), 3>& x, double t,                             \
      tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_SCALARS, (DataVector),
                        (CurvedScalarWave::Tags::Psi,
                         CurvedScalarWave::Tags::Pi))

#define INSTANTIATE_VECTORS(_, data)                                          \
  template tuples::TaggedTuple<TAG(data)> GaugeWaveConstantScalar::variables( \
      const tnsr::I<DTYPE(data), 3>& x, double t,                             \
      tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_VECTORS, (DataVector),
                        (CurvedScalarWave::Tags::Phi<3_st>))

#undef DTYPE
#undef TAG
#undef INSTANTIATE_SCALARS
#undef INSTANTIATE_VECTORS
}  // namespace ScalarTensor::Solutions
