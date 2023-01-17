// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Creators/CylindricalBinaryCompactObject.hpp"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include "Domain/BoundaryConditions/Periodic.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.tpp"
#include "Domain/CoordinateMaps/DiscreteRotation.hpp"
#include "Domain/CoordinateMaps/Interval.hpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/CoordinateMaps/ProductMaps.tpp"
#include "Domain/CoordinateMaps/TimeDependent/CubicScale.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ProductMaps.tpp"
#include "Domain/CoordinateMaps/TimeDependent/Rotation.hpp"
#include "Domain/CoordinateMaps/UniformCylindricalEndcap.hpp"
#include "Domain/CoordinateMaps/UniformCylindricalFlatEndcap.hpp"
#include "Domain/CoordinateMaps/UniformCylindricalSide.hpp"
#include "Domain/CoordinateMaps/Wedge.hpp"
#include "Domain/Creators/ExpandOverBlocks.hpp"
#include "Domain/DomainHelpers.hpp"
#include "Domain/FunctionsOfTime/FixedSpeedCubic.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "NumericalAlgorithms/RootFinding/QuadraticEquation.hpp"

namespace {
std::array<double, 3> rotate_to_z_axis(const std::array<double, 3> input) {
  return discrete_rotation(
      OrientationMap<3>{std::array<Direction<3>, 3>{Direction<3>::lower_zeta(),
                                                    Direction<3>::upper_eta(),
                                                    Direction<3>::upper_xi()}},
      input);
}
std::array<double, 3> flip_about_xy_plane(const std::array<double, 3> input) {
  return std::array<double, 3>{input[0], input[1], -input[2]};
}
}  // namespace

namespace domain::creators {
CylindricalBinaryCompactObject::CylindricalBinaryCompactObject(
    typename CenterA::type center_A, typename CenterB::type center_B,
    typename RadiusA::type radius_A, typename RadiusB::type radius_B,
    typename IncludeInnerSphereA::type include_inner_sphere_A,
    typename IncludeInnerSphereB::type include_inner_sphere_B,
    typename IncludeOuterSphere::type include_outer_sphere,
    typename OuterRadius::type outer_radius,
    const typename InitialRefinement::type& initial_refinement,
    const typename InitialGridPoints::type& initial_grid_points,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        inner_boundary_condition,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        outer_boundary_condition,
    const Options::Context& context)
    : center_A_(rotate_to_z_axis(center_A)),
      center_B_(rotate_to_z_axis(center_B)),
      radius_A_(radius_A),
      radius_B_(radius_B),
      include_inner_sphere_A_(include_inner_sphere_A),
      include_inner_sphere_B_(include_inner_sphere_B),
      include_outer_sphere_(include_outer_sphere),
      outer_radius_(outer_radius),
      inner_boundary_condition_(std::move(inner_boundary_condition)),
      outer_boundary_condition_(std::move(outer_boundary_condition)) {
  if (center_A_[2] <= 0.0) {
    PARSE_ERROR(
        context,
        "The x-coordinate of the input CenterA is expected to be positive");
  }
  if (center_B_[2] >= 0.0) {
    PARSE_ERROR(
        context,
        "The x-coordinate of the input CenterB is expected to be negative");
  }
  if (radius_A_ <= 0.0 or radius_B_ <= 0.0) {
    PARSE_ERROR(context, "RadiusA and RadiusB are expected to be positive");
  }
  if (radius_A_ < radius_B_) {
    PARSE_ERROR(context, "RadiusA should not be smaller than RadiusB");
  }
  if (std::abs(center_A_[2]) > std::abs(center_B_[2])) {
    PARSE_ERROR(context,
                "We expect |x_A| <= |x_B|, for x the x-coordinate of either "
                "CenterA or CenterB.  We should roughly have "
                "RadiusA x_A + RadiusB x_B = 0 (i.e. for BBHs the "
                "center of mass should be about at the origin).");
  }
  // The value 3.0 * (center_A_[2] - center_B_[2]) is what is
  // chosen in SpEC as the inner radius of the innermost outer sphere.
  if (outer_radius_ < 3.0 * (center_A_[2] - center_B_[2])) {
    PARSE_ERROR(context,
                "OuterRadius is too small. Please increase it "
                "beyond "
                    << 3.0 * (center_A_[2] - center_B_[2]));
  }

  if ((outer_boundary_condition_ == nullptr) xor
      (inner_boundary_condition_ == nullptr)) {
    PARSE_ERROR(context,
                "Must specify either both inner and outer boundary conditions "
                "or neither.");
  }
  using domain::BoundaryConditions::is_periodic;
  if (is_periodic(inner_boundary_condition_) or
      is_periodic(outer_boundary_condition_)) {
    PARSE_ERROR(
        context,
        "Cannot have periodic boundary conditions with a binary domain");
  }

  // The choices made below for the quantities xi, z_cutting_plane_,
  // and xi_min_sphere_e are the ones made in SpEC, and in the
  // Appendix of https://arxiv.org/abs/1206.3015.  Other choices could
  // be made that would still result in a reasonable Domain. In
  // particular, during a SpEC BBH evolution the excision boundaries
  // can sometimes get too close to z_cutting_plane_, and the
  // simulation must be halted and regridded with a different choice
  // of z_cutting_plane_, so it may be possible to choose a different
  // initial value of z_cutting_plane_ that reduces the number of such
  // regrids or eliminates them.

  // xi is the quantity in Eq. (A10) of
  // https://arxiv.org/abs/1206.3015 that represents how close the
  // cutting plane is to either center.  Unfortunately, there is a
  // discrepancy between what xi means in the paper and what it is in
  // the code.  I (Mark) think that this is a typo in the paper,
  // because otherwise the domain doesn't make sense.  To fix this,
  // either Eq. (A9) in the paper should have xi -> 1-xi, or Eq. (A10)
  // should have x_A and x_B swapped.
  // Here we will use the same definition of xi in Eq. (A10), but we
  // will swap xi -> 1-xi in Eq. (A9).
  // Therefore, xi = 0 means that the cutting plane passes through the center of
  // object B, and xi = 1 means that the cutting plane passes through
  // the center of object A.  Note that for |x_A| <= |x_B| (as assumed
  // above), xi is always <= 1/2.
  constexpr double xi_min = 0.25;
  // Same as Eq. (A10)
  const double xi =
      std::max(xi_min, std::abs(center_A_[2]) /
                           (std::abs(center_A_[2]) + std::abs(center_B_[2])));

  // Compute cutting plane
  // This is Eq. (A9) with xi -> 1-xi.
  z_cutting_plane_ = cut_spheres_offset_factor_ *
                     ((1.0 - xi) * center_B_[2] + xi * center_A_[2]);

  number_of_blocks_ = 46;
  if (include_inner_sphere_A) {
    number_of_blocks_ += 14;
  }
  if (include_inner_sphere_B) {
    number_of_blocks_ += 14;
  }
  if (include_outer_sphere) {
    number_of_blocks_ += 18;
  }

  // Add SphereE blocks if necessary.  Note that
  // https://arxiv.org/abs/1206.3015 has a mistake just above
  // Eq. (A.11) and the same mistake above Eq. (A.20), where it lists
  // the wrong mass ratio (for BBHs). The correct statement is that if
  // xi <= 1/3, this means that the mass ratio (for BBH) is large (>=2)
  // and we should add SphereE blocks.
  constexpr double xi_min_sphere_e = 1.0 / 3.0;
  if (xi <= xi_min_sphere_e) {
    // The following ERROR will be removed in an upcoming PR that
    // will support higher mass ratios.
    ERROR(
        "We currently only support domains where objects A and B are "
        "approximately the same size, and approximately the same distance from "
        "the origin.  More technically, we support xi > "
        << xi_min_sphere_e << ", but the value of xi is " << xi
        << ". Support for more general domains will be added in the near "
           "future");
  }

  // Create block names and groups
  auto add_filled_cylinder_name = [this](const std::string& prefix,
                                         const std::string& group_name) {
    for (const std::string& where :
         {"Center"s, "East"s, "North"s, "West"s, "South"s}) {
      const std::string name =
          std::string(prefix).append("FilledCylinder").append(where);
      block_names_.push_back(name);
      block_groups_[group_name].insert(name);
    }
  };
  auto add_cylinder_name = [this](const std::string& prefix,
                                  const std::string& group_name) {
    for (const std::string& where : {"East"s, "North"s, "West"s, "South"s}) {
      const std::string name =
          std::string(prefix).append("Cylinder").append(where);
      block_names_.push_back(name);
      block_groups_[group_name].insert(name);
    }
  };

  // CA Filled Cylinder
  // 5 blocks: 0 thru 4
  add_filled_cylinder_name("CA", "Outer");

  // CA Cylinder
  // 4 blocks: 5 thru 8
  add_cylinder_name("CA", "Outer");

  // EA Filled Cylinder
  // 5 blocks: 9 thru 13
  add_filled_cylinder_name("EA", "InnerA");

  // EA Cylinder
  // 4 blocks: 14 thru 17
  add_cylinder_name("EA", "InnerA");

  // EB Filled Cylinder
  // 5 blocks: 18 thru 22
  add_filled_cylinder_name("EB", "InnerB");

  // EB Cylinder
  // 4 blocks: 23 thru 26
  add_cylinder_name("EB", "InnerB");

  // MA Filled Cylinder
  // 5 blocks: 27 thru 31
  add_filled_cylinder_name("MA", "InnerA");

  // MB Filled Cylinder
  // 5 blocks: 32 thru 36
  add_filled_cylinder_name("MB", "InnerB");

  // CB Filled Cylinder
  // 5 blocks: 37 thru 41
  add_filled_cylinder_name("CB", "Outer");

  // CB Cylinder
  // 4 blocks: 42 thru 45
  add_cylinder_name("CB", "Outer");

  if (include_inner_sphere_A) {
    // 5 blocks
    add_filled_cylinder_name("InnerSphereEA", "InnerSphereA");
    // 5 blocks
    add_filled_cylinder_name("InnerSphereMA", "InnerSphereA");
    // 4 blocks
    add_cylinder_name("InnerSphereEA", "InnerSphereA");
  }
  if (include_inner_sphere_B) {
    // 5 blocks
    add_filled_cylinder_name("InnerSphereEB", "InnerSphereB");
    // 5 blocks
    add_filled_cylinder_name("InnerSphereMB", "InnerSphereB");
    // 4 blocks
    add_cylinder_name("InnerSphereEB", "InnerSphereB");
  }
  if (include_outer_sphere) {
    // 5 blocks
    add_filled_cylinder_name("OuterSphereCA", "OuterSphere");
    // 5 blocks
    add_filled_cylinder_name("OuterSphereCB", "OuterSphere");
    // 4 blocks
    add_cylinder_name("OuterSphereCA", "OuterSphere");
    // 4 blocks
    add_cylinder_name("OuterSphereCB", "OuterSphere");
  }

  // Expand initial refinement over all blocks
  const ExpandOverBlocks<size_t, 3> expand_over_blocks{block_names_,
                                                       block_groups_};
  try {
    initial_refinement_ = std::visit(expand_over_blocks, initial_refinement);
  } catch (const std::exception& error) {
    PARSE_ERROR(context, "Invalid 'InitialRefinement': " << error.what());
  }
  try {
    initial_grid_points_ = std::visit(expand_over_blocks, initial_grid_points);
  } catch (const std::exception& error) {
    PARSE_ERROR(context, "Invalid 'InitialGridPoints': " << error.what());
  }

  // Now we must change the initial refinement and initial grid points
  // for certain blocks, because the [r, theta, perp] directions do
  // not always correspond to [xi, eta, zeta].  The values in
  // initial_refinement_ must correspond to [xi, eta, zeta].
  //
  // In particular, for cylinders: [xi, eta, zeta] = [r, theta, perp]
  // but for filled cylinders: [xi, eta, zeta] = [perp, theta, r].

  auto swap_refinement_and_grid_points_xi_zeta = [this](const size_t block_id) {
    size_t val = gsl::at(initial_refinement_[block_id], 0);
    gsl::at(initial_refinement_[block_id], 0) =
        gsl::at(initial_refinement_[block_id], 2);
    gsl::at(initial_refinement_[block_id], 2) = val;
    val = gsl::at(initial_grid_points_[block_id], 0);
    gsl::at(initial_grid_points_[block_id], 0) =
        gsl::at(initial_grid_points_[block_id], 2);
    gsl::at(initial_grid_points_[block_id], 2) = val;
  };

  // CA Filled Cylinder
  // 5 blocks: 0 thru 4
  for (size_t block = 0; block < 5; ++block) {
    swap_refinement_and_grid_points_xi_zeta(block);
  }

  // EA Filled Cylinder
  // 5 blocks: 9 thru 13
  for (size_t block = 9; block < 14; ++block) {
    swap_refinement_and_grid_points_xi_zeta(block);
  }

  // EB Filled Cylinder
  // 5 blocks: 18 thru 22
  for (size_t block = 18; block < 23; ++block) {
    swap_refinement_and_grid_points_xi_zeta(block);
  }

  // MA Filled Cylinder
  // 5 blocks: 27 thru 31
  // MB Filled Cylinder
  // 5 blocks: 32 thru 36
  // CB Filled Cylinder
  // 5 blocks: 37 thru 41
  for (size_t block = 27; block < 42; ++block) {
    swap_refinement_and_grid_points_xi_zeta(block);
  }

  // Now do the filled cylinders for the inner and outer shells,
  // if they are present.
  size_t current_block = 46;
  if (include_inner_sphere_A) {
    for (size_t block = 0; block < 10; ++block) {
      swap_refinement_and_grid_points_xi_zeta(current_block++);
    }
    current_block += 4;
  }
  if (include_inner_sphere_B) {
    for (size_t block = 0; block < 10; ++block) {
      swap_refinement_and_grid_points_xi_zeta(current_block++);
    }
    current_block += 4;
  }
  if (include_outer_sphere) {
    for (size_t block = 0; block < 10; ++block) {
      swap_refinement_and_grid_points_xi_zeta(current_block++);
    }
  }
}

CylindricalBinaryCompactObject::CylindricalBinaryCompactObject(
    double initial_time, ExpansionMapOptions expansion_map_options,
    std::array<double, 3> initial_angular_velocity,
    typename CenterA::type center_A, typename CenterB::type center_B,
    typename RadiusA::type radius_A, typename RadiusB::type radius_B,
    typename IncludeInnerSphereA::type include_inner_sphere_A,
    typename IncludeInnerSphereB::type include_inner_sphere_B,
    typename IncludeOuterSphere::type include_outer_sphere,
    typename OuterRadius::type outer_radius,
    const typename InitialRefinement::type& initial_refinement,
    const typename InitialGridPoints::type& initial_grid_points,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        inner_boundary_condition,
    std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>
        outer_boundary_condition,
    const Options::Context& context)
    : CylindricalBinaryCompactObject(
          center_A, center_B, radius_A, radius_B, include_inner_sphere_A,
          include_inner_sphere_B, include_outer_sphere, outer_radius,
          initial_refinement, initial_grid_points,
          std::move(inner_boundary_condition),
          std::move(outer_boundary_condition), context) {
  is_time_dependent_ = true;
  initial_time_ = initial_time;
  expansion_map_options_ = expansion_map_options;
  initial_angular_velocity_ = initial_angular_velocity;
}

Domain<3> CylindricalBinaryCompactObject::create_domain() const {
  std::vector<std::unique_ptr<
      domain::CoordinateMapBase<Frame::BlockLogical, Frame::Inertial, 3>>>
      coordinate_maps{};

  const OrientationMap<3> rotate_to_x_axis{std::array<Direction<3>, 3>{
      Direction<3>::upper_zeta(), Direction<3>::upper_eta(),
      Direction<3>::lower_xi()}};

  const OrientationMap<3> rotate_to_minus_x_axis{std::array<Direction<3>, 3>{
      Direction<3>::lower_zeta(), Direction<3>::upper_eta(),
      Direction<3>::upper_xi()}};

  const std::array<double, 3> center_cutting_plane = {0.0, 0.0,
                                                      z_cutting_plane_};

  // The labels EA, EB, EE, etc are from Figure 20 of
  // https://arxiv.org/abs/1206.3015
  //
  // center_EA and radius_EA are the center and outer-radius of the
  // cylindered-sphere EA in Figure 20.
  //
  // center_EB and radius_EB are the center and outer-radius of the
  // cylindered-sphere EB in Figure 20.
  //
  // radius_MB is eq. A16 or A23 in the paper (depending on whether
  // the EE spheres exist), and is the radius of the circle where the EB
  // sphere intersects the cutting plane.
  const std::array<double, 3> center_EA = {
      0.0, 0.0, cut_spheres_offset_factor_ * center_A_[2]};
  const std::array<double, 3> center_EB = {
      0.0, 0.0, center_B_[2] * cut_spheres_offset_factor_};
  const double radius_MB =
      std::abs(cut_spheres_offset_factor_ * center_B_[2] - z_cutting_plane_);
  const double radius_EA =
      sqrt(square(center_EA[2] - z_cutting_plane_) + square(radius_MB));
  const double radius_EB =
      sqrt(2.0) * std::abs(center_EB[2] - z_cutting_plane_);

  // Construct vector<CoordMap>s that go from logical coordinates to
  // various blocks making up a unit right cylinder.  These blocks are
  // either the central square blocks, or the surrounding wedge
  // blocks. The radii and bounds are what are expected by the
  // UniformCylindricalEndcap maps, (except cylinder_inner_radius, which
  // determines the internal block boundaries inside the cylinder, and
  // which the UniformCylindricalEndcap maps don't care about).
  const double cylinder_inner_radius = 0.5;
  const double cylinder_outer_radius = 1.0;
  const double cylinder_lower_bound_z = -1.0;
  const double cylinder_upper_bound_z = 1.0;
  const auto logical_to_cylinder_center_maps =
      cyl_wedge_coord_map_center_blocks(cylinder_inner_radius,
                                        cylinder_lower_bound_z,
                                        cylinder_upper_bound_z, false);
  const auto logical_to_cylinder_surrounding_maps =
      cyl_wedge_coord_map_surrounding_blocks(
          cylinder_inner_radius, cylinder_outer_radius, cylinder_lower_bound_z,
          cylinder_upper_bound_z, false, 0.0);

  // Lambda that takes a UniformCylindricalEndcap map and a
  // DiscreteRotation map, composes it with the logical-to-cylinder
  // maps, and adds it to the list of coordinate maps. Also adds
  // boundary conditions if requested.
  auto add_endcap_to_list_of_maps =
      [&coordinate_maps, &logical_to_cylinder_center_maps,
       &logical_to_cylinder_surrounding_maps](
          const CoordinateMaps::UniformCylindricalEndcap& endcap_map,
          const CoordinateMaps::DiscreteRotation<3>& rotation_map) {
        auto new_logical_to_cylinder_center_maps =
            domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                    Frame::Inertial, 3>(
                logical_to_cylinder_center_maps, endcap_map, rotation_map);
        coordinate_maps.insert(
            coordinate_maps.end(),
            std::make_move_iterator(
                new_logical_to_cylinder_center_maps.begin()),
            std::make_move_iterator(new_logical_to_cylinder_center_maps.end()));
        auto new_logical_to_cylinder_surrounding_maps =
            domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                    Frame::Inertial, 3>(
                logical_to_cylinder_surrounding_maps, endcap_map, rotation_map);
        coordinate_maps.insert(
            coordinate_maps.end(),
            std::make_move_iterator(
                new_logical_to_cylinder_surrounding_maps.begin()),
            std::make_move_iterator(
                new_logical_to_cylinder_surrounding_maps.end()));
      };

  // Lambda that takes a UniformCylindricalFlatEndcap map and a
  // DiscreteRotation map, composes it with the logical-to-cylinder
  // maps, and adds it to the list of coordinate maps. Also adds
  // boundary conditions if requested.
  auto add_flat_endcap_to_list_of_maps =
      [&coordinate_maps, &logical_to_cylinder_center_maps,
       &logical_to_cylinder_surrounding_maps](
          const CoordinateMaps::UniformCylindricalFlatEndcap& endcap_map,
          const CoordinateMaps::DiscreteRotation<3>& rotation_map) {
        auto new_logical_to_cylinder_center_maps =
            domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                    Frame::Inertial, 3>(
                logical_to_cylinder_center_maps, endcap_map, rotation_map);
        coordinate_maps.insert(
            coordinate_maps.end(),
            std::make_move_iterator(
                new_logical_to_cylinder_center_maps.begin()),
            std::make_move_iterator(new_logical_to_cylinder_center_maps.end()));
        auto new_logical_to_cylinder_surrounding_maps =
            domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                    Frame::Inertial, 3>(
                logical_to_cylinder_surrounding_maps, endcap_map, rotation_map);
        coordinate_maps.insert(
            coordinate_maps.end(),
            std::make_move_iterator(
                new_logical_to_cylinder_surrounding_maps.begin()),
            std::make_move_iterator(
                new_logical_to_cylinder_surrounding_maps.end()));
      };

  // Construct vector<CoordMap>s that go from logical coordinates to
  // various blocks making up a right cylindrical shell of inner radius 1,
  // outer radius 2, and z-extents from -1 to +1.  These blocks are
  // either the central square blocks, or the surrounding wedge
  // blocks. The radii and bounds are what are expected by the
  // UniformCylindricalEndcap maps, (except cylinder_inner_radius, which
  // determines the internal block boundaries inside the cylinder, and
  // which the UniformCylindricalEndcap maps don't care about).
  const double cylindrical_shell_inner_radius = 1.0;
  const double cylindrical_shell_outer_radius = 2.0;
  const double cylindrical_shell_lower_bound_z = -1.0;
  const double cylindrical_shell_upper_bound_z = 1.0;
  const auto logical_to_cylindrical_shell_maps =
      cyl_wedge_coord_map_surrounding_blocks(
          cylindrical_shell_inner_radius, cylindrical_shell_outer_radius,
          cylindrical_shell_lower_bound_z, cylindrical_shell_upper_bound_z,
          false, 1.0);

  // Lambda that takes a UniformCylindricalSide map and a DiscreteRotation
  // map, composes it with the logical-to-cylinder maps, and adds it
  // to the list of coordinate maps.  Also adds boundary conditions if
  // requested.
  auto add_side_to_list_of_maps =
      [&coordinate_maps, &logical_to_cylindrical_shell_maps](
          const CoordinateMaps::UniformCylindricalSide& side_map,
          const CoordinateMaps::DiscreteRotation<3>& rotation_map) {
        auto new_logical_to_cylindrical_shell_maps =
            domain::make_vector_coordinate_map_base<Frame::BlockLogical,
                                                    Frame::Inertial, 3>(
                logical_to_cylindrical_shell_maps, side_map, rotation_map);
        coordinate_maps.insert(
            coordinate_maps.end(),
            std::make_move_iterator(
                new_logical_to_cylindrical_shell_maps.begin()),
            std::make_move_iterator(
                new_logical_to_cylindrical_shell_maps.end()));
      };

  // Inner radius of the outer C shell, if it exists.
  // If it doesn't exist, then it is the same as the outer_radius_.
  const double inner_radius_C = include_outer_sphere_
                                    ? 3.0 * (center_A_[2] - center_B_[2])
                                    : outer_radius_;

  // outer_radius_A is the outer radius of the inner sphere A, if it exists.
  // If the inner sphere A does not exist, then outer_radius_A is the same
  // as radius_A_.
  // If the inner sphere does exist, the algorithm for computing
  // outer_radius_A is the same as in SpEC when there is one inner shell.
  const double outer_radius_A =
      include_inner_sphere_A_
          ? radius_A_ +
                0.5 * (std::abs(z_cutting_plane_ - center_A_[2]) - radius_A_)
          : radius_A_;

  // outer_radius_B is the outer radius of the inner sphere B, if it exists.
  // If the inner sphere B does not exist, then outer_radius_B is the same
  // as radius_B_.
  // If the inner sphere does exist, the algorithm for computing
  // outer_radius_B is the same as in SpEC when there is one inner shell.
  const double outer_radius_B =
      include_inner_sphere_B_
          ? radius_B_ +
                0.5 * (std::abs(z_cutting_plane_ - center_B_[2]) - radius_B_)
          : radius_B_;

  // z_cut_CA_lower is the lower z_plane position for the CA endcap,
  // defined by https://arxiv.org/abs/1206.3015 in the bulleted list
  // after Eq. (A.19) EXCEPT that here we use a factor of 1.6 instead of 1.5
  // to put the plane farther from center_A.
  const double z_cut_CA_lower =
      z_cutting_plane_ + 1.6 * (center_EA[2] - z_cutting_plane_);
  // z_cut_CA_upper is the upper z_plane position for the CA endcap,
  // which isn't defined in https://arxiv.org/abs/1206.3015 (because the
  // maps are different).  We choose this plane to make the maps
  // less extreme.
  const double z_cut_CA_upper =
      std::max(0.5 * (z_cut_CA_lower + inner_radius_C), 0.7 * inner_radius_C);
  // z_cut_EA_upper is the upper z_plane position for the EA endcap,
  // which isn't defined in https://arxiv.org/abs/1206.3015 (because the
  // maps are different).  We choose this plane to make the maps
  // less extreme.
  const double z_cut_EA_upper = center_A_[2] + 0.7 * outer_radius_A;
  // z_cut_EA_lower is the lower z_plane position for the EA endcap,
  // which isn't defined in https://arxiv.org/abs/1206.3015 (because the
  // maps are different).  We choose this plane to make the maps
  // less extreme.
  const double z_cut_EA_lower = center_A_[2] - 0.7 * outer_radius_A;

  // CA Filled Cylinder
  // 5 blocks: 0 thru 4
  add_endcap_to_list_of_maps(
      CoordinateMaps::UniformCylindricalEndcap(center_EA, make_array<3>(0.0),
                                               radius_EA, inner_radius_C,
                                               z_cut_CA_lower, z_cut_CA_upper),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));

  // CA Cylinder
  // 4 blocks: 5 thru 8
  add_side_to_list_of_maps(
      CoordinateMaps::UniformCylindricalSide(
          // codecov complains about the next line being untested.
          // No idea why, since this entire function is called.
          // LCOV_EXCL_START
          center_EA, make_array<3>(0.0), radius_EA, inner_radius_C,
          // LCOV_EXCL_STOP
          z_cut_CA_lower, z_cutting_plane_, z_cut_CA_upper, z_cutting_plane_),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));

  // EA Filled Cylinder
  // 5 blocks: 9 thru 13
  add_endcap_to_list_of_maps(
      CoordinateMaps::UniformCylindricalEndcap(center_A_, center_EA,
                                               outer_radius_A, radius_EA,
                                               z_cut_EA_upper, z_cut_CA_lower),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));

  // EA Cylinder
  // 4 blocks: 14 thru 17
  add_side_to_list_of_maps(
      // For some reason codecov complains about the next line.
      CoordinateMaps::UniformCylindricalSide(  // LCOV_EXCL_LINE
          center_A_, center_EA, outer_radius_A, radius_EA, z_cut_EA_upper,
          z_cut_EA_lower, z_cut_CA_lower, z_cutting_plane_),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));

  // z_cut_CB_lower is the lower z_plane position for the CB endcap,
  // defined by https://arxiv.org/abs/1206.3015 in the bulleted list
  // after Eq. (A.19) EXCEPT that here we use a factor of 1.6 instead of 1.5
  // to put the plane farther from center_B.
  // Note here that 'lower' means 'farther from z=-infinity'
  // because we are on the -z side of the cutting plane.
  const double z_cut_CB_lower =
      z_cutting_plane_ + 1.6 * (center_EB[2] - z_cutting_plane_);
  // z_cut_CB_upper is the upper z_plane position for the CB endcap,
  // which isn't defined in https://arxiv.org/abs/1206.3015 (because the
  // maps are different).  We choose this plane to make the maps
  // less extreme. Note here that 'upper' means 'closer to z=-infinity'
  // because we are on the -z side of the cutting plane.
  const double z_cut_CB_upper =
      std::min(0.5 * (z_cut_CB_lower - inner_radius_C), -0.7 * inner_radius_C);
  // z_cut_EB_upper is the upper z_plane position for the EB endcap,
  // which isn't defined in https://arxiv.org/abs/1206.3015 (because the
  // maps are different).  We choose this plane to make the maps
  // less extreme.  Note here that 'upper' means 'closer to z=-infinity'
  // because we are on the -z side of the cutting plane.
  const double z_cut_EB_upper = center_B_[2] - 0.7 * outer_radius_B;
  // z_cut_EB_lower is the lower z_plane position for the EB endcap,
  // which isn't defined in https://arxiv.org/abs/1206.3015 (because the
  // maps are different).  We choose this plane to make the maps
  // less extreme. Note here that 'lower' means 'farther from z=-infinity'
  // because we are on the -z side of the cutting plane.
  const double z_cut_EB_lower = center_B_[2] + 0.7 * outer_radius_B;

  // EB Filled Cylinder
  // 5 blocks: 18 thru 22
  add_endcap_to_list_of_maps(
      CoordinateMaps::UniformCylindricalEndcap(
          flip_about_xy_plane(center_B_), flip_about_xy_plane(center_EB),
          outer_radius_B, radius_EB, -z_cut_EB_upper, -z_cut_CB_lower),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));

  // EB Cylinder
  // 4 blocks: 23 thru 26
  add_side_to_list_of_maps(
      CoordinateMaps::UniformCylindricalSide(
          flip_about_xy_plane(center_B_), flip_about_xy_plane(center_EB),
          outer_radius_B, radius_EB, -z_cut_EB_upper, -z_cut_EB_lower,
          -z_cut_CB_lower, -z_cutting_plane_),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));

  // MA Filled Cylinder
  // 5 blocks: 27 thru 31
  add_flat_endcap_to_list_of_maps(
      CoordinateMaps::UniformCylindricalFlatEndcap(
          flip_about_xy_plane(center_A_),
          flip_about_xy_plane(center_cutting_plane), outer_radius_A, radius_MB,
          -z_cut_EA_lower),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));
  // MB Filled Cylinder
  // 5 blocks: 32 thru 36
  add_flat_endcap_to_list_of_maps(
      // For some reason codecov complains about the next line.
      CoordinateMaps::UniformCylindricalFlatEndcap(  // LCOV_EXCL_LINE
          center_B_, center_cutting_plane, outer_radius_B, radius_MB,
          z_cut_EB_lower),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));

  // CB Filled Cylinder
  // 5 blocks: 37 thru 41
  add_endcap_to_list_of_maps(
      CoordinateMaps::UniformCylindricalEndcap(
          flip_about_xy_plane(center_EB), make_array<3>(0.0), radius_EB,
          inner_radius_C, -z_cut_CB_lower, -z_cut_CB_upper),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));

  // CB Cylinder
  // 4 blocks: 42 thru 45
  add_side_to_list_of_maps(
      CoordinateMaps::UniformCylindricalSide(
          flip_about_xy_plane(center_EB), make_array<3>(0.0), radius_EB,
          inner_radius_C, -z_cut_CB_lower, -z_cutting_plane_, -z_cut_CB_upper,
          -z_cutting_plane_),
      CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));

  if (include_inner_sphere_A_) {
    const double z_cut_upper = center_A_[2] + 0.7 * radius_A_;
    const double z_cut_lower = center_A_[2] - 0.7 * radius_A_;
    // InnerSphereEA Filled Cylinder
    // 5 blocks
    add_endcap_to_list_of_maps(
        // For some reason codecov complains about the next function.
        // LCOV_EXCL_START
        CoordinateMaps::UniformCylindricalEndcap(center_A_, center_A_,
                                                 radius_A_, outer_radius_A,
                                                 z_cut_upper, z_cut_EA_upper),
        // LCOV_EXCL_START
        CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));
    // InnerSphereMA Filled Cylinder
    // 5 blocks
    add_endcap_to_list_of_maps(
        CoordinateMaps::UniformCylindricalEndcap(
            flip_about_xy_plane(center_A_), flip_about_xy_plane(center_A_),
            radius_A_, outer_radius_A, -z_cut_lower, -z_cut_EA_lower),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));
    // InnerSphereEA Cylinder
    // 4 blocks
    add_side_to_list_of_maps(
        // For some reason codecov complains about the next line.
        CoordinateMaps::UniformCylindricalSide(  // LCOV_EXCL_LINE
            center_A_, center_A_, radius_A_, outer_radius_A, z_cut_upper,
            z_cut_lower, z_cut_EA_upper, z_cut_EA_lower),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));
  }
  if (include_inner_sphere_B_) {
    // Note here that 'upper' means 'closer to z=-infinity'
    // because we are on the -z side of the cutting plane.
    const double z_cut_upper = center_B_[2] - 0.7 * radius_B_;
    const double z_cut_lower = center_B_[2] + 0.7 * radius_B_;
    // InnerSphereEB Filled Cylinder
    // 5 blocks
    add_endcap_to_list_of_maps(
        CoordinateMaps::UniformCylindricalEndcap(
            flip_about_xy_plane(center_B_), flip_about_xy_plane(center_B_),
            radius_B_, outer_radius_B, -z_cut_upper, -z_cut_EB_upper),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));
    // InnerSphereMB Filled Cylinder
    // 5 blocks
    add_endcap_to_list_of_maps(
        // For some reason codecov complains about the next function.
        // LCOV_EXCL_START
        CoordinateMaps::UniformCylindricalEndcap(center_B_, center_B_,
                                                 radius_B_, outer_radius_B,
                                                 z_cut_lower, z_cut_EB_lower),
        // LCOV_EXCL_STOP
        CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));
    // InnerSphereEB Cylinder
    // 4 blocks
    add_side_to_list_of_maps(
        CoordinateMaps::UniformCylindricalSide(
            flip_about_xy_plane(center_B_), flip_about_xy_plane(center_B_),
            radius_B_, outer_radius_B, -z_cut_upper, -z_cut_lower,
            -z_cut_EB_upper, -z_cut_EB_lower),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));
  }
  if (include_outer_sphere_) {
    const double z_cut_CA_outer = 0.7 * outer_radius_;
    const double z_cut_CB_outer = -0.7 * outer_radius_;
    // OuterCA Filled Cylinder
    // 5 blocks
    add_endcap_to_list_of_maps(
        CoordinateMaps::UniformCylindricalEndcap(
            make_array<3>(0.0), make_array<3>(0.0), inner_radius_C,
            outer_radius_, z_cut_CA_upper, z_cut_CA_outer),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));
    // OuterCB Filled Cylinder
    // 5 blocks
    add_endcap_to_list_of_maps(
        CoordinateMaps::UniformCylindricalEndcap(
            make_array<3>(0.0), make_array<3>(0.0), inner_radius_C,
            outer_radius_, -z_cut_CB_upper, -z_cut_CB_outer),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));
    // OuterCA Cylinder
    // 4 blocks
    add_side_to_list_of_maps(
        CoordinateMaps::UniformCylindricalSide(
            make_array<3>(0.0), make_array<3>(0.0), inner_radius_C,
            outer_radius_, z_cut_CA_upper, z_cutting_plane_, z_cut_CA_outer,
            z_cutting_plane_),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_x_axis));
    // OuterCB Cylinder
    // 4 blocks
    add_side_to_list_of_maps(
        CoordinateMaps::UniformCylindricalSide(
            make_array<3>(0.0), make_array<3>(0.0), inner_radius_C,
            outer_radius_, -z_cut_CB_upper, -z_cutting_plane_, -z_cut_CB_outer,
            -z_cutting_plane_),
        CoordinateMaps::DiscreteRotation<3>(rotate_to_minus_x_axis));
  }

  Domain<3> domain{std::move(coordinate_maps)};

  if (is_time_dependent_) {
    // Same map for all the blocks for now.
    // When size, shape, cutx, skew are added then we will have different
    // maps for some blocks.
    using CubicScaleMap = domain::CoordinateMaps::TimeDependent::CubicScale<3>;
    using CubicScaleMapForComposition =
        domain::CoordinateMap<Frame::Grid, Frame::Inertial, CubicScaleMap>;

    using RotationMap3D = domain::CoordinateMaps::TimeDependent::Rotation<3>;
    using RotationMapForComposition =
        domain::CoordinateMap<Frame::Grid, Frame::Inertial, RotationMap3D>;

    using CubicScaleAndRotationMapForComposition =
        domain::CoordinateMap<Frame::Grid, Frame::Inertial, CubicScaleMap,
                              RotationMap3D>;
    std::unique_ptr<domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>>
        full_map = std::make_unique<CubicScaleAndRotationMapForComposition>(
            domain::push_back(
                CubicScaleMapForComposition{CubicScaleMap{
                    outer_radius_, expansion_function_of_time_name_,
                    expansion_function_of_time_name_ + "OuterBoundary"s}},
                RotationMapForComposition{
                    RotationMap3D{rotation_function_of_time_name_}}));

    for (size_t block = 0; block < number_of_blocks_; ++block) {
      domain.inject_time_dependent_map_for_block(block, full_map->get_clone());
    }
  }
  return domain;
}

std::vector<DirectionMap<
    3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
CylindricalBinaryCompactObject::external_boundary_conditions() const {
  if (outer_boundary_condition_ == nullptr) {
    return {};
  }
  std::vector<DirectionMap<
      3, std::unique_ptr<domain::BoundaryConditions::BoundaryCondition>>>
      boundary_conditions{number_of_blocks_};
  for (size_t i = 0; i < 5; ++i) {
    if (not include_outer_sphere_) {
      // CA Filled Cylinder
      boundary_conditions[i][Direction<3>::upper_zeta()] =
          outer_boundary_condition_->get_clone();
      // CB Filled Cylinder
      boundary_conditions[i + 37][Direction<3>::upper_zeta()] =
          outer_boundary_condition_->get_clone();
    }
    if (not include_inner_sphere_A_) {
      // EA Filled Cylinder
      boundary_conditions[i + 9][Direction<3>::lower_zeta()] =
          inner_boundary_condition_->get_clone();
      // MA Filled Cylinder
      boundary_conditions[i + 27][Direction<3>::lower_zeta()] =
          inner_boundary_condition_->get_clone();
    }
    if (not include_inner_sphere_B_) {
      // EB Filled Cylinder
      boundary_conditions[i + 18][Direction<3>::lower_zeta()] =
          inner_boundary_condition_->get_clone();
      // MB Filled Cylinder
      boundary_conditions[i + 32][Direction<3>::lower_zeta()] =
          inner_boundary_condition_->get_clone();
    }
  }
  for (size_t i = 0; i < 4; ++i) {
    if (not include_outer_sphere_) {
      // CA Cylinder
      boundary_conditions[i + 5][Direction<3>::upper_xi()] =
          outer_boundary_condition_->get_clone();
      // CB Cylinder
      boundary_conditions[i + 42][Direction<3>::upper_xi()] =
          outer_boundary_condition_->get_clone();
    }
    if (not include_inner_sphere_A_) {
      // EA Cylinder
      boundary_conditions[i + 14][Direction<3>::lower_xi()] =
          inner_boundary_condition_->get_clone();
    }
    if (not include_inner_sphere_B_) {
      // EB Cylinder
      boundary_conditions[i + 23][Direction<3>::lower_xi()] =
          inner_boundary_condition_->get_clone();
    }
  }

  size_t last_block = 46;
  if (include_inner_sphere_A_) {
    for (size_t i = 0; i < 5; ++i) {
      // InnerSphereEA Filled Cylinder
      boundary_conditions[last_block + i][Direction<3>::lower_zeta()] =
          inner_boundary_condition_->get_clone();
      // InnerSphereMA Filled Cylinder
      boundary_conditions[last_block + i + 5][Direction<3>::lower_zeta()] =
          inner_boundary_condition_->get_clone();
    }
    for (size_t i = 0; i < 4; ++i) {
      // InnerSphereEA Cylinder
      boundary_conditions[last_block + i + 10][Direction<3>::lower_xi()] =
          inner_boundary_condition_->get_clone();
    }
    last_block += 14;
  }
  if (include_inner_sphere_B_) {
    for (size_t i = 0; i < 5; ++i) {
      // InnerSphereEB Filled Cylinder
      boundary_conditions[last_block + i][Direction<3>::lower_zeta()] =
          inner_boundary_condition_->get_clone();
      // InnerSphereMB Filled Cylinder
      boundary_conditions[last_block + i + 5][Direction<3>::lower_zeta()] =
          inner_boundary_condition_->get_clone();
    }
    for (size_t i = 0; i < 4; ++i) {
      // InnerSphereEB Cylinder
      boundary_conditions[last_block + i + 10][Direction<3>::lower_xi()] =
          inner_boundary_condition_->get_clone();
    }
    last_block += 14;
  }
  if (include_outer_sphere_) {
    for (size_t i = 0; i < 5; ++i) {
      // OuterCA Filled Cylinder
      boundary_conditions[last_block + i][Direction<3>::upper_zeta()] =
          outer_boundary_condition_->get_clone();
      // OuterCB Filled Cylinder
      boundary_conditions[last_block + i + 5][Direction<3>::upper_zeta()] =
          outer_boundary_condition_->get_clone();
    }
    for (size_t i = 0; i < 4; ++i) {
      // OuterCA Cylinder
      boundary_conditions[last_block + i + 10][Direction<3>::upper_xi()] =
          outer_boundary_condition_->get_clone();
      // OuterCB Cylinder
      boundary_conditions[last_block + i + 14][Direction<3>::upper_xi()] =
          outer_boundary_condition_->get_clone();
    }
  }
  return boundary_conditions;
}

std::vector<std::array<size_t, 3>>
CylindricalBinaryCompactObject::initial_extents() const {
  return initial_grid_points_;
}

std::vector<std::array<size_t, 3>>
CylindricalBinaryCompactObject::initial_refinement_levels() const {
  return initial_refinement_;
}

std::unordered_map<std::string,
                   std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
CylindricalBinaryCompactObject::functions_of_time(
    const std::unordered_map<std::string, double>& initial_expiration_times)
    const {
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      result{};
  if (not is_time_dependent_) {
    return result;
  }

  // Get existing function of time names that are used for the maps and assign
  // their initial expiration time to infinity (i.e. not expiring)
  std::unordered_map<std::string, double> expiration_times{
      {expansion_function_of_time_name_,
       std::numeric_limits<double>::infinity()},
      {rotation_function_of_time_name_,
       std::numeric_limits<double>::infinity()}};

  // If we have control systems, overwrite these expiration times with the ones
  // supplied by the control system
  for (auto& [name, expr_time] : initial_expiration_times) {
    expiration_times[name] = expr_time;
  }

  // ExpansionMap FunctionOfTime for the function \f$a(t)\f$ in the
  // domain::CoordinateMaps::TimeDependent::CubicScale map
  result[expansion_function_of_time_name_] =
      std::make_unique<FunctionsOfTime::PiecewisePolynomial<2>>(
          initial_time_,
          std::array<DataVector, 3>{{{expansion_map_options_.initial_data[0]},
                                     {expansion_map_options_.initial_data[1]},
                                     {0.0}}},
          expiration_times.at(expansion_function_of_time_name_));

  // ExpansionMap FunctionOfTime for the function \f$b(t)\f$ in the
  // domain::CoordinateMaps::TimeDependent::CubicScale map
  result[expansion_function_of_time_name_ + "OuterBoundary"s] =
      std::make_unique<FunctionsOfTime::FixedSpeedCubic>(
          // codecov seems to not recognize that the next line is actually
          // run during the tests.
          // LCOV_EXCL_START
          1.0, initial_time_, expansion_map_options_.outer_boundary_velocity,
          // LCOV_EXCL_STOP
          expansion_map_options_.outer_boundary_decay_time);

  // RotationMap FunctionOfTime for the rotation angles about each
  // axis.  The initial rotation angles don't matter as we never
  // actually use the angles themselves. We only use their derivatives
  // (omega) to determine map parameters. In theory we could determine
  // each initital angle from the input axis-angle representation, but
  // we don't need to.
  result[rotation_function_of_time_name_] =
      std::make_unique<FunctionsOfTime::QuaternionFunctionOfTime<3>>(
          initial_time_,
          std::array<DataVector, 1>{DataVector{1.0, 0.0, 0.0, 0.0}},
          std::array<DataVector, 4>{
              {{3, 0.0},
               {initial_angular_velocity_[0], initial_angular_velocity_[1],
                initial_angular_velocity_[2]},
               {3, 0.0},
               {3, 0.0}}},
          expiration_times.at(rotation_function_of_time_name_));
  return result;
}
}  // namespace domain::creators
