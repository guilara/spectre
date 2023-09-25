// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Protocols/Mutator.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Initialization/InitialData.hpp"
#include "Evolution/Initialization/Tags.hpp"
#include "Evolution/Systems/CurvedScalarWave/BackgroundSpacetime.hpp"
#include "Evolution/Systems/CurvedScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/TagsDeclarations.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "PointwiseFunctions/AnalyticData/Tags.hpp"
#include "Utilities/ProtocolHelpers.hpp"

namespace fe::ScalarDriver::Initialization {
/// \ingroup InitializationGroup
/// \brief Mutator meant to be used with
/// `Initialization::Actions::AddSimpleTags` to initialize the constraint
/// damping parameters of the CurvedScalarWave system
///
/// DataBox changes:
/// - Adds:
///   * `fe::ScalarDriver::Tags::ConstraintGamma1`
///   * `fe::ScalarDriver::Tags::ConstraintGamma2`
/// - Removes: nothing
/// - Modifies: nothing

template <size_t Dim>
struct InitializeConstraintDampingGammas
    : tt::ConformsTo<db::protocols::Mutator> {
  using return_tags =
      tmpl::list<Tags::ConstraintGamma1, Tags::ConstraintGamma2>;
  using argument_tags = tmpl::list<domain::Tags::Mesh<Dim>>;

  static void apply(const gsl::not_null<Scalar<DataVector>*> gamma1,
                    const gsl::not_null<Scalar<DataVector>*> gamma2,
                    const Mesh<Dim>& mesh) {
    const size_t number_of_grid_points = mesh.number_of_grid_points();
    *gamma1 = Scalar<DataVector>{number_of_grid_points, 0.};
    *gamma2 = Scalar<DataVector>{number_of_grid_points, 1.};
  }
};

}  // namespace fe::ScalarDriver::Initialization
