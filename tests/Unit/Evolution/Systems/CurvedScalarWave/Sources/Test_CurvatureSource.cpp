// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <string>
#include <type_traits>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Determinant.hpp"
#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/CoordinateMaps/Affine.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/ProductMaps.tpp"
#include "Evolution/Systems/CurvedScalarWave/Sources/CurvatureSource.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/PointwiseFunctions/AnalyticSolutions/GeneralRelativity/VerifyGrSolution.hpp"
#include "Helpers/PointwiseFunctions/AnalyticSolutions/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.tpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/GeneralRelativity/Ricci.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace {
using Affine = domain::CoordinateMaps::Affine;
using Affine3D = domain::CoordinateMaps::ProductOf3Maps<Affine, Affine, Affine>;

template <typename FrameType, typename DataType>
tnsr::I<DataType, 3, FrameType> spatial_coords(const DataType& used_for_size) {
  auto x = make_with_value<tnsr::I<DataType, 3, FrameType>>(used_for_size, 0.0);
  get<0>(x) = 1.32;
  get<1>(x) = 0.82;
  get<2>(x) = 1.24;
  return x;
}

template <typename FrameType, typename DataType>
void test_compute_scalar_curvature_source(const DataType& used_for_size) {
  // We want to compare with analytic results for the Schwarzschild metric
  // The Kretchmann for a non-rotating black hole in both Schwarzschild and
  // Kerr-Schild coordinates is 48 M^2 / r^6 , where r^2 = x^2 + y^2 + z^2.

  // Add the following libraries to the RunSingleTest cmake file:
  //   LinearOperators
  //   CurvedScalarWave
  //   CurvedScalarWaveSources
  //   GeneralRelativitySolutions

  // Define solution parameters
  const double mass = 1.0;
  const std::array<double, 3> spin{{0.0, 0.0, 0.0}};
  const std::array<double, 3> center{{0.0, 0.0, 0.0}};
  // Create instance of solution
  const gr::Solutions::KerrSchild& solution{mass, spin, center};

  // Setup grid
  const size_t num_points_1d = 8;
  const std::array<double, 3> lower_bound{{0.8, 1.22, 1.30}};
  const std::array<double, 3> upper_bound{{0.82, 1.24, 1.32}};
  const size_t SpatialDim = 3;
  Mesh<SpatialDim> mesh{num_points_1d, Spectral::Basis::Legendre,
                        Spectral::Quadrature::GaussLobatto};
  const auto coord_map =
      domain::make_coordinate_map<Frame::ElementLogical, FrameType>(Affine3D{
          Affine{-1., 1., lower_bound[0], upper_bound[0]},
          Affine{-1., 1., lower_bound[1], upper_bound[1]},
          Affine{-1., 1., lower_bound[2], upper_bound[2]},
      });
  const size_t num_points_3d = num_points_1d * num_points_1d * num_points_1d;
  // Setup coordinates
  const auto x_logical = logical_coordinates(mesh);
  const auto x = coord_map(x_logical);
  // Arbitrary time for time-independent solution.
  const double t = std::numeric_limits<double>::signaling_NaN();

  const auto vars = solution.variables(
      x, t, typename gr::Solutions::KerrSchild::tags<DataType, FrameType>{});

  const auto& g = get<gr::Tags::SpatialMetric<3, FrameType, DataType>>(vars);
  const auto& ig =
      get<gr::Tags::InverseSpatialMetric<3, FrameType, DataType>>(vars);
  const auto& sqrt_det_spatial_metric =
      get<typename gr::Tags::SqrtDetSpatialMetric<DataType>>(vars);
  const auto& extrinsic_curvature =
      get<typename gr::Tags::ExtrinsicCurvature<3, FrameType, DataType>>(vars);
  const auto& spatial_christoffel_second_kind = get<
      typename gr::Tags::SpatialChristoffelSecondKind<3, FrameType, DataType>>(
      vars);

  // Compute needed partial derivatives
  const auto& deriv_spatial_christoffel_second_kind = partial_derivative(
      spatial_christoffel_second_kind, mesh, coord_map.inv_jacobian(x_logical));
  const auto& deriv_extrinsic_curvature = partial_derivative(
      extrinsic_curvature, mesh, coord_map.inv_jacobian(x_logical));

  // Compute spatial Ricci tensor
  const auto& spatial_ricci_tensor =
      gr::ricci_tensor<3, FrameType, IndexType::Spatial, DataType>(
          spatial_christoffel_second_kind,
          deriv_spatial_christoffel_second_kind);
  // Compute spatial Ricci scalar
  const auto& spatial_ricci_scalar =
      gr::ricci_scalar<3, FrameType, IndexType::Spatial, DataType>(
          spatial_ricci_tensor, ig);

  // Compute Weyl Electric tensor
  const auto& weyl_electric_tensor = gr::weyl_electric<3, FrameType, DataType>(
      spatial_ricci_tensor, extrinsic_curvature, ig);
  // Compute Weyl Electric Square scalar
  const auto& weyl_electric_scalar =
      gr::weyl_electric_scalar<3, FrameType, DataType>(weyl_electric_tensor,
                                                       ig);

  // Compute Weyl Magnetic tensor
  const auto& weyl_magnetic_tensor = gr::weyl_magnetic<FrameType, DataType>(
      deriv_extrinsic_curvature, g, sqrt_det_spatial_metric);
  // Compute Weyl Magnetic Square scalar
  const auto& weyl_magnetic_scalar =
      gr::weyl_magnetic_scalar<FrameType, DataType>(weyl_magnetic_tensor, ig);

  // Compute source term
  const double first_coupling_psi = 1.0;
  const double second_coupling_psi = 0.0;
  const double mass_psi = 0.0;
  const auto& psi = make_with_value<Scalar<DataVector>>(num_points_3d, 0.0);

  const auto& source_term =
      CurvedScalarWave::Sources::compute_scalar_curvature_source(
          weyl_electric_scalar, weyl_magnetic_scalar, psi, first_coupling_psi,
          second_coupling_psi, mass_psi);
}

SPECTRE_TEST_CASE(
    "Unit.Evolution.Systems.CurvedScalarWave.Sources.CurvatureSource",
    "[Unit][Evolution][Options]") {
  test_compute_scalar_curvature_source<::Frame::Inertial, DataVector>(
      DataVector(5));
}

}  // namespace
