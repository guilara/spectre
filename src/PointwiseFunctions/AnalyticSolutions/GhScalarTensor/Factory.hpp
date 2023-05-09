// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/GaugeWaveConstantScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/GaugeWaveScalarWave.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/KerrSchildScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/KerrSphericalHarmonic.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/MinkowskiZeroScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/MinkowskiScalarWave.hpp"
#include "Utilities/TMPL.hpp"

namespace GeneralizedHarmonic::Solutions {
/// ScalarTensor solutions wrapped for GH
namespace ScalarTensor {
/// \brief List of all analytic solutions
using all_solutions =
    tmpl::list<GeneralizedHarmonic::Solutions::WrappedGr<
                   ::ScalarTensor::Solutions::MinkowskiZeroScalar>,
               GeneralizedHarmonic::Solutions::WrappedGr<
                   ::ScalarTensor::Solutions::KerrSchildScalar>,
               GeneralizedHarmonic::Solutions::WrappedGr<
                   ::ScalarTensor::Solutions::GaugeWaveConstantScalar>,
               GeneralizedHarmonic::Solutions::WrappedGr<
                   ::ScalarTensor::Solutions::MinkowskiScalarWave>,
               GeneralizedHarmonic::Solutions::WrappedGr<
                   ::ScalarTensor::Solutions::GaugeWaveScalarWave>,
               GeneralizedHarmonic::Solutions::WrappedGr<
                   ::ScalarTensor::Solutions::KerrSphericalHarmonic>>;
}  // namespace ScalarTensor
}  // namespace GeneralizedHarmonic::Solutions
