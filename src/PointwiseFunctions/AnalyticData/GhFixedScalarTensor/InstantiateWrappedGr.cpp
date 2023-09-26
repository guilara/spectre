// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/DataBox/Prefixes.hpp"
#include "PointwiseFunctions/AnalyticData/FixedScalarTensor/FixedDecoupledScalar/KerrSphericalHarmonic.hpp"
#include "PointwiseFunctions/AnalyticData/FixedScalarTensor/FixedDecoupledScalar/MinkowskiZeroScalar.hpp"
#include "PointwiseFunctions/AnalyticData/FixedScalarTensor/FixedDecoupledScalar/SphericalKerrSchildSH.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.tpp"
#include "Utilities/GenerateInstantiations.hpp"

GENERATE_INSTANTIATIONS(
    WRAPPED_GR_INSTANTIATE,
    (fe::DecoupledScalar::AnalyticData::MinkowskiZeroScalar,
     fe::DecoupledScalar::AnalyticData::KerrSphericalHarmonic,
     fe::DecoupledScalar::AnalyticData::SphericalKerrSchildSH))
