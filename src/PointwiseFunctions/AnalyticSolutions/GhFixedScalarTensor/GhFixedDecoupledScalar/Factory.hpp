// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticSolutions/FixedScalarTensor/FixedDecoupledScalar/MinkowskiZeroScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/FixedScalarTensor/FixedDecoupledScalar/KerrSchildScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/FixedScalarTensor/FixedDecoupledScalar/KerrSphericalHarmonic.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Utilities/TMPL.hpp"

namespace gh::Solutions {
/// ScalarTensor solutions wrapped for GH
namespace fe::DecoupledScalar {
/// \brief List of all analytic solutions
using all_solutions = tmpl::list<gh::Solutions::WrappedGr<
    ::fe::DecoupledScalar::Solutions::MinkowskiZeroScalar>,
    gh::Solutions::WrappedGr<
    ::fe::DecoupledScalar::Solutions::KerrSchildScalar>,
    gh::Solutions::WrappedGr<
    ::fe::DecoupledScalar::Solutions::KerrSphericalHarmonic>>;
}  // namespace fe::DecoupledScalar
}  // namespace gh::Solutions
