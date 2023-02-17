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

#include <sstream>
#include "Parallel/Printf.hpp"
#include "Utilities/GetOutput.hpp"

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
  const double mass = 0.5;
  const std::array<double, 3> spin{{0.1, 0.2, 0.3}};
  const std::array<double, 3> center{{1.0, 3.0, 2.0}};
  // Create instance of wrapped solution
  const GeneralizedHarmonic::Solutions::WrappedGr<gr::Solutions::KerrSchild>&
      wrapped_ks_solution{mass, spin, center};
  const gr::Solutions::KerrSchild& ks_solution{mass, spin, center};

  const DataVector data_vector(10, 1.0);
  Parallel::printf("KerrSchild::volume_dim = %d \n",
                   gr::Solutions::KerrSchild::volume_dim);
  const tnsr::I<DataVector, gr::Solutions::KerrSchild::volume_dim,
                Frame::Inertial>
      x{data_vector};
  const double t = 1.0;

  const auto wrapped_gh_vars = wrapped_ks_solution.variables(
      x, t,
      tmpl::list<
          gr::Tags::SpacetimeMetric<gr::Solutions::KerrSchild::volume_dim,
                                    Frame::Inertial, DataVector>,
          GeneralizedHarmonic::Tags::Pi<gr::Solutions::KerrSchild::volume_dim,
                                        Frame::Inertial>,
          GeneralizedHarmonic::Tags::Phi<gr::Solutions::KerrSchild::volume_dim,
                                         Frame::Inertial>>{});

  const auto ks_vars = ks_solution.variables(
      x, t,
      tmpl::list<gr::Tags::SpatialChristoffelFirstKind<3, ::Frame::Inertial,
                                                       DataVector>,
                 gr::Tags::SpatialChristoffelSecondKind<3, ::Frame::Inertial,
                                                        DataVector>,
                 gr::Tags::TraceSpatialChristoffelSecondKind<
                     3, ::Frame::Inertial, DataVector>>{});

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

  const auto& christ_down_trace =
      get<gr::Tags::TraceSpatialChristoffelSecondKind<3, ::Frame::Inertial,
                                                      DataVector>>(ks_vars);
  // const DataVector christ_data = christ_down_trace.get();

  Parallel::printf("Size of spacetime_metric: %d \n", spacetime_metric.size());
  Parallel::printf("Size of pi : %d \n", pi.size());
  Parallel::printf("Size of phi : %d \n", phi.size());
  Parallel::printf("Size of christ_down_trace : %d \n",
                   christ_down_trace.size());
  // for (size_t index = 0; index < spacetime_metric.size(); ++index) {
  //   Parallel::printf("index: %d \n", index);
  //   Parallel::printf("spacetime_metric[index]: %lf \n",
  //                    spacetime_metric[index]);
  //   Parallel::printf("pi[index]: %lf \n", pi[index]);
  // }
  // for (size_t index = 0; index < christ_down_trace.size(); ++index) {
  //   Parallel::printf("christ_down_trace[index]: %lf \n",
  //                    christ_down_trace[index]);
  // }

  // const DataVector check_vec{0.0, 0.0, 0.0};
  const auto check_tensor =
      make_with_value<tnsr::I<DataVector, 3>>(christ_down_trace.size(), 0.0);
  const auto check_tensor_2 =
      make_with_value<tnsr::aa<DataVector, 3, ::Frame::Inertial>>(
          data_vector.size(), 0.0);
  // Parallel::printf("check size: %d \n", check_vec.size());

  std::ostringstream os;

  os << get_output(check_tensor_2);

  Parallel::printf(os.str());

  for (size_t index = 0; index < 3; ++index) {
  }

  // CHECK_ITERABLE_APPROX(christ_down_trace, check_tensor);
  // CHECK_ITERABLE_APPROX(spacetime_metric, check_tensor_2);

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
