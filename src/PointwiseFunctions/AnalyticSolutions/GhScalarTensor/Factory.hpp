// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/KerrSchildScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/KerrSphericalHarmonic.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/MinkowskiZeroScalar.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Solutions {
/// ScalarTensor solutions wrapped for GH
namespace ScalarTensor {
/// \brief List of all analytic solutions
using all_solutions = tmpl::list<
    gh::Solutions::WrappedGr<::ScalarTensor::Solutions::MinkowskiZeroScalar>,
    gh::Solutions::WrappedGr<::ScalarTensor::Solutions::KerrSchildScalar>,
    gh::Solutions::WrappedGr<::ScalarTensor::Solutions::KerrSphericalHarmonic>>;
}  // namespace ScalarTensor
}  // namespace gh::Solutions
