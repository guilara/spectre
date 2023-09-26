// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "PointwiseFunctions/AnalyticData/FixedScalarTensor/FixedDecoupledScalar/KerrSphericalHarmonic.hpp"
#include "PointwiseFunctions/AnalyticData/FixedScalarTensor/FixedDecoupledScalar/MinkowskiZeroScalar.hpp"
#include "PointwiseFunctions/AnalyticData/FixedScalarTensor/FixedDecoupledScalar/SphericalKerrSchildSH.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Utilities/TMPL.hpp"

/// FixedScalarTensor solutions wrapped for GH
namespace gh::fe::DecoupledScalar::AnalyticData {

/// \brief List of all analytic solutions
using all_analytic_data =
    tmpl::list<gh::Solutions::WrappedGr<
                   ::fe::DecoupledScalar::AnalyticData::MinkowskiZeroScalar>,
               gh::Solutions::WrappedGr<
                   ::fe::DecoupledScalar::AnalyticData::KerrSphericalHarmonic>,
               gh::Solutions::WrappedGr<
                   ::fe::DecoupledScalar::AnalyticData::SphericalKerrSchildSH>>;

}  // namespace gh::fe::DecoupledScalar::AnalyticData
