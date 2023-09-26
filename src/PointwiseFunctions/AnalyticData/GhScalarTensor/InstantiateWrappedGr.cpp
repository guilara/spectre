// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/DataBox/Prefixes.hpp"
#include "PointwiseFunctions/AnalyticData/ScalarTensor/KerrSphericalHarmonic.hpp"
#include "PointwiseFunctions/AnalyticData/ScalarTensor/MinkowskiZeroScalar.hpp"
#include "PointwiseFunctions/AnalyticData/ScalarTensor/SphericalKerrSchildSH.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.tpp"
#include "Utilities/GenerateInstantiations.hpp"

GENERATE_INSTANTIATIONS(WRAPPED_GR_INSTANTIATE,
                        (ScalarTensor::AnalyticData::MinkowskiZeroScalar,
                         ScalarTensor::AnalyticData::KerrSphericalHarmonic,
                         ScalarTensor::AnalyticData::SphericalKerrSchildSH))
