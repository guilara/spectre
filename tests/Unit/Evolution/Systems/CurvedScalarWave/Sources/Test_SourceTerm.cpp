// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <limits>
#include <string>

#include "DataStructures/DataBox/Prefixes.hpp"  // IWYU pragma: keep
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/BackgroundSpacetime.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/GaugeWave.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/SphericalKerrSchild.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/WrappedGr.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Phi.hpp"
#include "PointwiseFunctions/GeneralRelativity/GeneralizedHarmonic/Pi.hpp"
#include "PointwiseFunctions/GeneralRelativity/SpacetimeMetric.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace {

void test_compute_scalar_source () {
//
}

// We want to test WrappedGr
// Add the following libraries to the RunSingleTest cmake file:
//   CurvedScalarWave
//   CurvedScalarWaveSources
//   GeneralRelativitySolutions
//
void test_background__spacetime () {
    // Define solution parameters
    const double mass = 1.0;
    const std::array<double, 3> spin{{0.0, 0.0, 0.0}};
    const std::array<double, 3> center{{0.0, 0.0, 0.0}};
    // Create instance of wrapped solution
    const GeneralizedHarmonic::Solutions::WrappedGr<gr::Solutions::KerrSchild>&
        wrapped_ks_solution{mass, spin, center};

}

SPECTRE_TEST_CASE("Unit.Evolution.Systems.CurvedScalarWave.Sources.SourceTerm",
                  "[Unit][Evolution][Options]") {
    test_compute_scalar_source();
    test_background__spacetime();
}

} // namespace
