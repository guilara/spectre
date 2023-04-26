// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/DataBox/Prefixes.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.tpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/KerrSchildScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/MinkowskiZeroScalar.hpp"
#include "Utilities/GenerateInstantiations.hpp"

GENERATE_INSTANTIATIONS(WRAPPED_GR_INSTANTIATE,
                        (ScalarTensor::Solutions::MinkowskiZeroScalar,
                         ScalarTensor::Solutions::KerrSchildScalar))
