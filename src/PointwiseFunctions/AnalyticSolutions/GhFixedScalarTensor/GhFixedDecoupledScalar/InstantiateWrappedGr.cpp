// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/DataBox/Prefixes.hpp"
#include "PointwiseFunctions/AnalyticSolutions/FixedScalarTensor/FixedDecoupledScalar/MinkowskiZeroScalar.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "Utilities/GenerateInstantiations.hpp"

GENERATE_INSTANTIATIONS(WRAPPED_GR_INSTANTIATE,
                        (fe::DecoupledScalar::Solutions::MinkowskiZeroScalar))
