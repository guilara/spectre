// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Helpers/PointwiseFunctions/AnalyticSolutions/GeneralRelativity/CheckWrappedGrConsistency.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/GaugeWaveConstantScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/KerrSchildScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/MinkowskiZeroScalar.hpp"

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.GhScalarTensor.WrappedGr",
    "[Unit][PointwiseFunctions]") {
  const double time = 2.0;
  const tnsr::I<DataVector, 3, Frame::Inertial> coords{DataVector{3.0, 4.0}};

  check_wrapped_gr_solution_consistency(
      GeneralizedHarmonic::Solutions::WrappedGr<
          ScalarTensor::Solutions::MinkowskiZeroScalar>{1.0},
          ScalarTensor::Solutions::MinkowskiZeroScalar{1.0}, coords, time);

  check_wrapped_gr_solution_consistency(
      GeneralizedHarmonic::Solutions::WrappedGr<
          ScalarTensor::Solutions::KerrSchildScalar>{1.0, 2.0},
      ScalarTensor::Solutions::KerrSchildScalar{1.0, 2.0}, coords, time);

  check_wrapped_gr_solution_consistency(
      GeneralizedHarmonic::Solutions::WrappedGr<
          ScalarTensor::Solutions::GaugeWaveConstantScalar>{0.5, 1.0, 1.0},
      ScalarTensor::Solutions::GaugeWaveConstantScalar{0.5, 1.0, 1.0}, coords,
      time);
}
