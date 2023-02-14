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

#include "Parallel/Printf.hpp"

namespace {

void test_spacetime_riemann_tensor() {
  //
}

// We want to test WrappedGr
// Add the following libraries to the RunSingleTest cmake file:
//   CurvedScalarWave
//   CurvedScalarWaveSources
//   GeneralRelativitySolutions
//
void test_background__spacetime() {
  // Define solution parameters
  const double mass = 1.0;
  const std::array<double, 3> spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> center{{0.0, 0.0, 0.0}};
  // Create instance of wrapped solution
  const GeneralizedHarmonic::Solutions::WrappedGr<gr::Solutions::KerrSchild>&
      wrapped_ks_solution{mass, spin, center};

  const DataVector data_vector{3.0, 4.0};
  const tnsr::I<DataVector, gr::Solutions::KerrSchild::volume_dim,
                Frame::Inertial>
      x{data_vector};
  const double t = 44.44;

  const auto wrapped_gh_vars = wrapped_ks_solution.variables(
      x, t,
      tmpl::list<
          gr::Tags::SpacetimeMetric<gr::Solutions::KerrSchild::volume_dim,
                                    Frame::Inertial, DataVector>,
          GeneralizedHarmonic::Tags::Pi<gr::Solutions::KerrSchild::volume_dim,
                                        Frame::Inertial>,
          GeneralizedHarmonic::Tags::Phi<gr::Solutions::KerrSchild::volume_dim,
                                         Frame::Inertial>>{});

  const auto& spacetime_metric =
      get<gr::Tags::SpacetimeMetric<gr::Solutions::KerrSchild::volume_dim,
                                    Frame::Inertial, DataVector>>(
          wrapped_gh_vars);
  const auto& pi =
      get<GeneralizedHarmonic::Tags::Pi<gr::Solutions::KerrSchild::volume_dim,
                                        Frame::Inertial>>(wrapped_gh_vars);
  const auto& phi =
      get<GeneralizedHarmonic::Tags::Phi<gr::Solutions::KerrSchild::volume_dim,
                                        Frame::Inertial>>(wrapped_gh_vars);

  Parallel::printf("Size of spacetime_metric: %d \n", spacetime_metric.size());
  Parallel::printf("Size of pi : %d \n", pi.size());
  Parallel::printf("Size of phi : %d \n", phi.size());

  //   auto test_riemann_tensor =
  //       spacetime_riemann_tensor<DataVector, 4_st, ::Frame::Inertial>(
  //           pi, phi, deriv_pi, deriv_phi, inverse_spatial_metric);
}

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.GeneralRelativity.SpacetimeRiemann.",
                 "[PointwiseFunctions][Unit]") {
  test_spacetime_riemann_tensor();
  test_background__spacetime();
}

}  // namespace
