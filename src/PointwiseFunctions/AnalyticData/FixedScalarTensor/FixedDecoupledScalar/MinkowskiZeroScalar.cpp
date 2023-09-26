// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/FixedScalarTensor/FixedDecoupledScalar/MinkowskiZeroScalar.hpp"

#include <cmath>
#include <cstddef>

#include "DataStructures/DataVector.hpp"     // IWYU pragma: keep
#include "DataStructures/Tensor/Tensor.hpp"  // IWYU pragma: keep
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

// IWYU pragma:  no_include "DataStructures/Tensor/TypeAliases.hpp"

namespace fe::DecoupledScalar::AnalyticData {

MinkowskiZeroScalar::MinkowskiZeroScalar(double amplitude,
                                         double amplitude_driver) {
  amplitude_ = amplitude;
  amplitude_driver_ = amplitude_driver;
}

std::unique_ptr<evolution::initial_data::InitialData>
MinkowskiZeroScalar::get_clone() const {
  return std::make_unique<MinkowskiZeroScalar>(*this);
}

MinkowskiZeroScalar::MinkowskiZeroScalar(CkMigrateMessage* msg)
    : InitialData(msg) {}

void MinkowskiZeroScalar::pup(PUP::er& p) {
  InitialData::pup(p);
  p | amplitude_;
  p | amplitude_driver_;
  p | background_spacetime_;
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Psi> MinkowskiZeroScalar::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Psi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, amplitude_)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Phi<3_st>>
MinkowskiZeroScalar::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Phi<3_st>> /*meta*/) const {
  return {make_with_value<tnsr::i<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<CurvedScalarWave::Tags::Pi> MinkowskiZeroScalar::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<CurvedScalarWave::Tags::Pi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, 0.0)};
}

//

template <typename DataType>
tuples::TaggedTuple<fe::ScalarDriver::Tags::Psi> MinkowskiZeroScalar::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<fe::ScalarDriver::Tags::Psi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, amplitude_driver_)};
}

template <typename DataType>
tuples::TaggedTuple<fe::ScalarDriver::Tags::Phi<3_st>>
MinkowskiZeroScalar::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<fe::ScalarDriver::Tags::Phi<3_st>> /*meta*/) const {
  return {make_with_value<tnsr::i<DataType, 3>>(x, 0.0)};
}

template <typename DataType>
tuples::TaggedTuple<fe::ScalarDriver::Tags::Pi> MinkowskiZeroScalar::variables(
    const tnsr::I<DataType, 3>& x,
    tmpl::list<fe::ScalarDriver::Tags::Pi> /*meta*/) const {
  return {make_with_value<Scalar<DataType>>(x, 0.0)};
}

PUP::able::PUP_ID MinkowskiZeroScalar::my_PUP_ID = 0;

bool operator==(const MinkowskiZeroScalar& lhs,
                const MinkowskiZeroScalar& rhs) {
  return lhs.background_spacetime_ == rhs.background_spacetime_ and
         lhs.amplitude_ == rhs.amplitude_ and
         lhs.amplitude_driver_ == rhs.amplitude_driver_;
}

bool operator!=(const MinkowskiZeroScalar& lhs,
                const MinkowskiZeroScalar& rhs) {
  return not(lhs == rhs);
}

#define DTYPE(data) BOOST_PP_TUPLE_ELEM(0, data)
#define TAG(data) BOOST_PP_TUPLE_ELEM(1, data)

#define INSTANTIATE_SCALARS(_, data)                                      \
  template tuples::TaggedTuple<TAG(data)> MinkowskiZeroScalar::variables( \
      const tnsr::I<DTYPE(data), 3>& x, tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_SCALARS, (DataVector),
                        (CurvedScalarWave::Tags::Psi,
                         CurvedScalarWave::Tags::Pi,
                         fe::ScalarDriver::Tags::Psi,
                         fe::ScalarDriver::Tags::Pi))

#define INSTANTIATE_VECTORS(_, data)                                      \
  template tuples::TaggedTuple<TAG(data)> MinkowskiZeroScalar::variables( \
      const tnsr::I<DTYPE(data), 3>& x, tmpl::list<TAG(data)> /*meta*/) const;

GENERATE_INSTANTIATIONS(INSTANTIATE_VECTORS, (DataVector),
                        (CurvedScalarWave::Tags::Phi<3_st>,
                         fe::ScalarDriver::Tags::Phi<3_st>))

#undef DTYPE
#undef TAG
#undef INSTANTIATE_SCALARS
#undef INSTANTIATE_VECTORS
}  // namespace fe::DecoupledScalar::AnalyticData
