// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticData/ScalarTensor/KerrSphericalHarmonic.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::AnalyticData {
/// ScalarTensor solutions wrapped for GH
namespace ScalarTensor {
/// \brief List of all analytic solutions
using all_analytic_data = tmpl::list<gh::Solutions::WrappedGr<
    ::ScalarTensor::AnalyticData::KerrSphericalHarmonic>>;
}  // namespace ScalarTensor
}  // namespace gh::AnalyticData
