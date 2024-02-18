// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/SphereTimeDependentMaps.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <variant>

#include "DataStructures/DataVector.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/SphereTransition.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/FunctionsOfTime/SettleToConstant.hpp"
#include "Domain/FunctionsOfTime/SettleToConstantQuaternion.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrHorizon.hpp"
#include "Utilities/ErrorHandling/Error.hpp"

namespace domain::creators::sphere {

TimeDependentMapOptions::TimeDependentMapOptions(
    const double initial_time, const ShapeMapOptions& shape_map_options,
    const RotationMapOptions& rotation_map_options,
    const ExpansionMapOptions& expansion_map_options,
    const TranslationMapOptions& translation_map_options)
    : initial_time_(initial_time),
      shape_map_options_(shape_map_options),
      rotation_map_options_(rotation_map_options),
      expansion_map_options_(expansion_map_options),
      translation_map_options_(translation_map_options) {}

std::unordered_map<std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
TimeDependentMapOptions::create_functions_of_time(
    const double inner_radius,
    const std::unordered_map<std::string, double>& initial_expiration_times)
    const {
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      result{};

  // Get existing function of time names that are used for the maps and assign
  // their initial expiration time to infinity (i.e. not expiring)
  std::unordered_map<std::string, double> expiration_times{
      {size_name, std::numeric_limits<double>::infinity()},
      {shape_name, std::numeric_limits<double>::infinity()},
      {expansion_name, std::numeric_limits<double>::infinity()},
      {expansion_outer_name, std::numeric_limits<double>::infinity()},
      {rotation_name, std::numeric_limits<double>::infinity()},
      {translation_name, std::numeric_limits<double>::infinity()}};

  // If we have control systems, overwrite these expiration times with the ones
  // supplied by the control system
  for (const auto& [name, expr_time] : initial_expiration_times) {
    expiration_times[name] = expr_time;
  }

  DataVector shape_zeros{
      ylm::Spherepack::spectral_size(shape_map_options_.l_max,
                                     shape_map_options_.l_max),
      0.0};
  DataVector shape_func{};
  DataVector size_func{1, 0.0};

  if (shape_map_options_.initial_values.has_value()) {
    if (std::holds_alternative<KerrSchildFromBoyerLindquist>(
            shape_map_options_.initial_values.value())) {
      const ylm::Spherepack ylm{shape_map_options_.l_max,
                                shape_map_options_.l_max};
      const auto& mass_and_spin = std::get<KerrSchildFromBoyerLindquist>(
          shape_map_options_.initial_values.value());
      const DataVector radial_distortion =
          1.0 - get(gr::Solutions::kerr_schild_radius_from_boyer_lindquist(
                    inner_radius, ylm.theta_phi_points(), mass_and_spin.mass,
                    mass_and_spin.spin)) /
                    inner_radius;
      shape_func = ylm.phys_to_spec(radial_distortion);
      // Transform from SPHEREPACK to actual Ylm for size func
      size_func[0] = shape_func[0] * sqrt(0.5 * M_PI);
      // Set l=0 for shape map to 0 because size is going to be used
      shape_func[0] = 0.0;
    }
  } else {
    shape_func = shape_zeros;
    size_func[0] = 0.0;
  }

  // ShapeMap FunctionOfTime
  result[shape_name] =
      std::make_unique<FunctionsOfTime::PiecewisePolynomial<2>>(
          initial_time_,
          std::array<DataVector, 3>{
              {std::move(shape_func), shape_zeros, shape_zeros}},
          expiration_times.at(shape_name));

  DataVector size_deriv{1, 0.0};
  DataVector size_2nd_deriv{1, 0.0};

  // Size FunctionOfTime (used in ShapeMap)
  result[size_name] = std::make_unique<FunctionsOfTime::PiecewisePolynomial<3>>(
      initial_time_,
      std::array<DataVector, 4>{{std::move(size_func),
                                 std::move(size_deriv),
                                 std::move(size_2nd_deriv),
                                 {0.0}}},
      expiration_times.at(size_name));

  // ExpansionMap FunctionOfTime
  // Note: Missing expiration times in SettleToConstant
  result[expansion_name] = std::make_unique<FunctionsOfTime::SettleToConstant>(
      std::array<DataVector, 3>{
          {{gsl::at(expansion_map_options_.initial_values, 0)},
           {gsl::at(expansion_map_options_.initial_values, 1)},
           {0.0}}},
      initial_time_, expansion_map_options_.decay_timescale);

  // ExpansionMap in the Outer regionFunctionOfTime
  result[expansion_outer_name] =
      std::make_unique<FunctionsOfTime::PiecewisePolynomial<2>>(
          initial_time_, std::array<DataVector, 3>{{{1.0}, {0.0}, {0.0}}},
          expiration_times.at(expansion_outer_name));

  // RotationMap FunctionOfTime
  // Note: Missing expiration times in SettleToConstantQuaternion
  result[rotation_name] =
      std::make_unique<FunctionsOfTime::SettleToConstantQuaternion>(
          std::array<DataVector, 3>{
              DataVector{1.0, 0.0, 0.0, 0.0},
              DataVector{
                  0.0,
                  gsl::at(rotation_map_options_.initial_angular_velocity, 0),
                  gsl::at(rotation_map_options_.initial_angular_velocity, 1),
                  gsl::at(rotation_map_options_.initial_angular_velocity, 2)},
              {4, 0.0}},
          initial_time_, rotation_map_options_.decay_timescale);

  DataVector initial_translation_center_temp{3, 0.0};
  DataVector initial_translation_velocity_temp{3, 0.0};
  for (size_t i = 0; i < 3; i++) {
    initial_translation_center_temp[i] =
        gsl::at(translation_map_options_.initial_values.front(), i);
    initial_translation_velocity_temp[i] =
        gsl::at(translation_map_options_.initial_values.back(), i);
  }

  // Translation FunctionOfTime
  result[translation_name] =
      std::make_unique<FunctionsOfTime::PiecewisePolynomial<2>>(
          initial_time_,
          std::array<DataVector, 3>{
              {std::move(initial_translation_center_temp),
               std::move(initial_translation_velocity_temp),
               {3, 0.0}}},
          expiration_times.at(translation_name));

  return result;
}

void TimeDependentMapOptions::build_maps(
    const std::array<double, 3>& center, const double inner_radius,
    const double outer_radius,
    std::pair<double, double> translation_transition_radii) {
  std::unique_ptr<domain::CoordinateMaps::ShapeMapTransitionFunctions::
                      ShapeMapTransitionFunction>
      transition_func =
          std::make_unique<domain::CoordinateMaps::ShapeMapTransitionFunctions::
                               SphereTransition>(inner_radius, outer_radius);
  shape_map_ = ShapeMap{center,
                        shape_map_options_.l_max,
                        shape_map_options_.l_max,
                        std::move(transition_func),
                        shape_name,
                        size_name};

  inner_rot_scale_trans_map_ = RotScaleTransMap{
      std::pair<std::string, std::string>{expansion_name, expansion_outer_name},
      rotation_name,
      translation_name,
      translation_transition_radii.first,
      translation_transition_radii.second,
      domain::CoordinateMaps::TimeDependent::RotScaleTrans<
          3>::BlockRegion::Inner};

  transition_rot_scale_trans_map_ = RotScaleTransMap{
      std::pair<std::string, std::string>{expansion_name, expansion_outer_name},
      rotation_name,
      translation_name,
      translation_transition_radii.first,
      translation_transition_radii.second,
      domain::CoordinateMaps::TimeDependent::RotScaleTrans<
          3>::BlockRegion::Transition};
}

// If you edit any of the functions below, be sure to update the documentation
// in the Sphere domain creator as well as this class' documentation.
TimeDependentMapOptions::MapType<Frame::Distorted, Frame::Inertial>
TimeDependentMapOptions::distorted_to_inertial_map(
    const bool include_distorted_map) const {
  if (include_distorted_map) {
    return std::make_unique<DistortedToInertialComposition>(
        inner_rot_scale_trans_map_);
  } else {
    return nullptr;
  }
}

TimeDependentMapOptions::MapType<Frame::Grid, Frame::Distorted>
TimeDependentMapOptions::grid_to_distorted_map(
    const bool include_distorted_map) const {
  if (include_distorted_map) {
    return std::make_unique<GridToDistortedComposition>(shape_map_);
  } else {
    return nullptr;
  }
}

TimeDependentMapOptions::MapType<Frame::Grid, Frame::Inertial>
TimeDependentMapOptions::grid_to_inertial_map(
    const bool include_distorted_map, const bool use_rigid_translation) const {
  if (include_distorted_map) {
    return std::make_unique<GridToInertialComposition>(
        shape_map_, inner_rot_scale_trans_map_);
  } else if (use_rigid_translation) {
    return std::make_unique<GridToInertialSimple>(inner_rot_scale_trans_map_);
  } else {
    return std::make_unique<GridToInertialSimple>(
        transition_rot_scale_trans_map_);
  }
}
}  // namespace domain::creators::sphere
