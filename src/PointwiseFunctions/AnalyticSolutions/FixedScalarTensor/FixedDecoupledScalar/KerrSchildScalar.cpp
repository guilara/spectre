// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/FixedScalarTensor/FixedDecoupledScalar/KerrSchildScalar.hpp"

#include <array>
#include <cmath>
#include <cstddef>

#include "DataStructures/DataVector.hpp"     // IWYU pragma: keep
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedDecoupledScalar/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

// IWYU pragma:  no_include "DataStructures/Tensor/TypeAliases.hpp"

namespace fe::DecoupledScalar::Solutions {

KerrSchildScalar::KerrSchildScalar(double mass, double amplitude) {
  mass_ = mass;
  amplitude_ = amplitude;
  background_spacetime_ =
      gr::Solutions::KerrSchild{// BH mass
                                mass_,
                                // Dimensionless spin
                                std::array<double, 3>{{0.0, 0.0, 0.0}},
                                // Center
                                std::array<double, 3>{{0.0, 0.0, 0.0}}};
}

std::unique_ptr<evolution::initial_data::InitialData>
KerrSchildScalar::get_clone() const {
  return std::make_unique<KerrSchildScalar>(*this);
}

KerrSchildScalar::KerrSchildScalar(CkMigrateMessage* msg) : InitialData(msg) {}

void KerrSchildScalar::pup(PUP::er& p) {
  InitialData::pup(p);
  p | mass_;
  p | amplitude_;
  p | background_spacetime_;
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Psi> KerrSchildScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<CurvedScalarWave::Tags::Psi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, amplitude_)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Phi<3_st>>
KerrSchildScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<CurvedScalarWave::Tags::Phi<3_st>> /*meta*/) const {
  return {make_with_value<tnsr::i<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Pi> KerrSchildScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<CurvedScalarWave::Tags::Pi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, 0.0)};
}

//

template <typename DataType>
tuples::TaggedTuple<fe::ScalarDriver::Tags::Psi> KerrSchildScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<fe::ScalarDriver::Tags::Psi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, amplitude_)};
}

template <typename DataType>
tuples::TaggedTuple<fe::ScalarDriver::Tags::Phi<3_st>>
KerrSchildScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<fe::ScalarDriver::Tags::Phi<3_st>> /*meta*/) const {
  return {make_with_value<tnsr::i<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<fe::ScalarDriver::Tags::Pi> KerrSchildScalar::variables(
    const tnsr::I<DataType, 3>& x, double /*t*/,
    tmpl::list<fe::ScalarDriver::Tags::Pi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, 0.0)};
}

PUP::able::PUP_ID KerrSchildScalar::my_PUP_ID = 0;

bool operator==(const KerrSchildScalar& lhs, const KerrSchildScalar& rhs) {
  return lhs.background_spacetime_ == rhs.background_spacetime_ and
         lhs.mass_ == rhs.mass_ and lhs.amplitude_ == rhs.amplitude_;
}

bool operator!=(const KerrSchildScalar& lhs, const KerrSchildScalar& rhs) {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE_SCALARS(_, data)                                   \
  template tuples::TaggedTuple<TAG(data)> KerrSchildScalar::variables( \
      const tnsr::I<DTYPE(data), 3>& x, double t,                      \
      tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_SCALARS, (DataVector),
                        (CurvedScalarWave::Tags::Psi,
                         CurvedScalarWave::Tags::Pi,
                         fe::ScalarDriver::Tags::Psi,
                         fe::ScalarDriver::Tags::Pi))

#define INSTANTIATE_VECTORS(_, data)                                   \
  template tuples::TaggedTuple<TAG(data)> KerrSchildScalar::variables( \
      const tnsr::I<DTYPE(data), 3>& x, double t,                      \
      tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_VECTORS, (DataVector),
                        (CurvedScalarWave::Tags::Phi<3_st>,
                         fe::ScalarDriver::Tags::Phi<3_st>))

#undef DTYPE
#undef TAG
#undef INSTANTIATE_SCALARS
#undef INSTANTIATE_VECTORS
}  // namespace fe::DecoupledScalar::Solutions
