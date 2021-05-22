// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <type_traits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Expressions/Evaluate.hpp"
#include "DataStructures/Tensor/Expressions/TensorExpression.hpp"
#include "DataStructures/Tensor/IndexType.hpp"
#include "DataStructures/Tensor/Symmetry.hpp"
#include "Helpers/DataStructures/Tensor/Expressions/EvaluateRank0TestHelpers.hpp"
#include "Helpers/DataStructures/Tensor/Expressions/EvaluateRank1TestHelpers.hpp"
#include "Helpers/DataStructures/Tensor/Expressions/EvaluateRank2TestHelpers.hpp"
#include "Helpers/DataStructures/Tensor/Expressions/EvaluateRank3TestHelpers.hpp"
#include "Helpers/DataStructures/Tensor/Expressions/EvaluateRank4TestHelpers.hpp"

namespace {
template <auto&... TensorIndices>
void test_contains_indices_to_contract_impl(const bool expected) {
  CHECK(TensorExpressions::detail::contains_indices_to_contract<sizeof...(
            TensorIndices)>(
            {{std::decay_t<decltype(TensorIndices)>::value...}}) == expected);
}

void test_contains_indices_to_contract() {
  test_contains_indices_to_contract_impl<ti_a, ti_b, ti_c>(false);
  test_contains_indices_to_contract_impl<ti_I, ti_j>(false);
  test_contains_indices_to_contract_impl<ti_j>(false);
  test_contains_indices_to_contract_impl(false);

  test_contains_indices_to_contract_impl<ti_d, ti_D>(true);
  test_contains_indices_to_contract_impl<ti_I, ti_i>(true);
  test_contains_indices_to_contract_impl<ti_a, ti_K, ti_B, ti_b>(true);
  test_contains_indices_to_contract_impl<ti_j, ti_c, ti_J, ti_A, ti_a>(true);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.DataStructures.Tensor.Expression.Evaluate",
                  "[DataStructures][Unit]") {
  test_contains_indices_to_contract();

  // Rank 0: double
  TestHelpers::TensorExpressions::test_evaluate_rank_0<double>(-7.31);

  // Rank 0: DataVector
  TestHelpers::TensorExpressions::test_evaluate_rank_0<DataVector>(
      DataVector{-3.1, 9.4, 0.0, -3.1, 2.4, 9.8});

  // Rank 1: double; spacetime
  TestHelpers::TensorExpressions::test_evaluate_rank_1<double, SpacetimeIndex,
                                                       UpLo::Lo, ti_a>();
  TestHelpers::TensorExpressions::test_evaluate_rank_1<double, SpacetimeIndex,
                                                       UpLo::Lo, ti_b>();
  TestHelpers::TensorExpressions::test_evaluate_rank_1<double, SpacetimeIndex,
                                                       UpLo::Up, ti_A>();
  TestHelpers::TensorExpressions::test_evaluate_rank_1<double, SpacetimeIndex,
                                                       UpLo::Up, ti_B>();

  // Rank 1: double; spatial
  TestHelpers::TensorExpressions::test_evaluate_rank_1<double, SpatialIndex,
                                                       UpLo::Lo, ti_i>();
  TestHelpers::TensorExpressions::test_evaluate_rank_1<double, SpatialIndex,
                                                       UpLo::Lo, ti_j>();
  TestHelpers::TensorExpressions::test_evaluate_rank_1<double, SpatialIndex,
                                                       UpLo::Up, ti_I>();
  TestHelpers::TensorExpressions::test_evaluate_rank_1<double, SpatialIndex,
                                                       UpLo::Up, ti_J>();

  // Rank 1: DataVector
  TestHelpers::TensorExpressions::test_evaluate_rank_1<DataVector, SpatialIndex,
                                                       UpLo::Up, ti_L>();

  // Rank 2: double; nonsymmetric; spacetime only
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Lo, UpLo::Lo, ti_a, ti_b>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Up, UpLo::Up, ti_A, ti_B>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Lo, UpLo::Lo, ti_d, ti_c>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Up, UpLo::Up, ti_D, ti_C>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Lo, UpLo::Up, ti_e, ti_F>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Up, UpLo::Lo, ti_F, ti_e>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Lo, UpLo::Up, ti_g, ti_B>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Up, UpLo::Lo, ti_G, ti_b>();

  // Rank 2: double; nonsymmetric; spatial only
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpatialIndex, UpLo::Lo, UpLo::Lo, ti_i, ti_j>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpatialIndex, UpLo::Up, UpLo::Up, ti_I, ti_J>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpatialIndex, UpLo::Lo, UpLo::Lo, ti_j, ti_i>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpatialIndex, UpLo::Up, UpLo::Up, ti_J, ti_I>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpatialIndex, UpLo::Lo, UpLo::Up, ti_i, ti_J>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpatialIndex, UpLo::Up, UpLo::Lo, ti_I, ti_j>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpatialIndex, UpLo::Lo, UpLo::Up, ti_j, ti_I>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpatialIndex, UpLo::Up, UpLo::Lo, ti_J, ti_i>();

  // Rank 2: double; nonsymmetric; spacetime and spatial mixed
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpatialIndex, UpLo::Lo, UpLo::Up, ti_c, ti_I>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpatialIndex, UpLo::Up, UpLo::Lo, ti_A, ti_i>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpacetimeIndex, UpLo::Up, UpLo::Lo, ti_J, ti_a>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpacetimeIndex, UpLo::Lo, UpLo::Up, ti_i, ti_B>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpatialIndex, UpLo::Lo, UpLo::Lo, ti_e, ti_j>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpacetimeIndex, UpLo::Lo, UpLo::Lo, ti_i, ti_d>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpacetimeIndex, SpatialIndex, UpLo::Up, UpLo::Up, ti_C, ti_I>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      double, SpatialIndex, SpacetimeIndex, UpLo::Up, UpLo::Up, ti_J, ti_A>();

  // Rank 2: double; symmetric; spacetime
  TestHelpers::TensorExpressions::test_evaluate_rank_2_symmetric<
      double, SpacetimeIndex, UpLo::Lo, ti_a, ti_d>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_symmetric<
      double, SpacetimeIndex, UpLo::Up, ti_G, ti_B>();

  // Rank 2: double; symmetric; spatial
  TestHelpers::TensorExpressions::test_evaluate_rank_2_symmetric<
      double, SpatialIndex, UpLo::Lo, ti_j, ti_i>();
  TestHelpers::TensorExpressions::test_evaluate_rank_2_symmetric<
      double, SpatialIndex, UpLo::Up, ti_I, ti_J>();

  // Rank 2: DataVector; nonsymmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_2_no_symmetry<
      DataVector, SpacetimeIndex, SpacetimeIndex, UpLo::Lo, UpLo::Up, ti_f,
      ti_G>();

  // Rank 2: DataVector; symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_2_symmetric<
      DataVector, SpatialIndex, UpLo::Lo, ti_j, ti_i>();

  // Rank 3: double; nonsymmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_no_symmetry<
      double, SpacetimeIndex, SpatialIndex, SpacetimeIndex, UpLo::Up, UpLo::Lo,
      UpLo::Up, ti_D, ti_j, ti_B>();

  // Rank 3: double; first and second indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_ab_symmetry<
      double, SpacetimeIndex, SpacetimeIndex, UpLo::Lo, UpLo::Up, ti_b, ti_a,
      ti_C>();

  // Rank 3: double; first and third indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_ac_symmetry<
      double, SpatialIndex, SpacetimeIndex, UpLo::Lo, UpLo::Lo, ti_i, ti_f,
      ti_j>();

  // Rank 3: double; second and third indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_bc_symmetry<
      double, SpacetimeIndex, SpatialIndex, UpLo::Lo, UpLo::Up, ti_d, ti_J,
      ti_I>();

  // Rank 3: double; symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_abc_symmetry<
      double, SpacetimeIndex, UpLo::Lo, ti_f, ti_d, ti_a>();

  // Rank 3: DataVector; nonsymmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_no_symmetry<
      DataVector, SpacetimeIndex, SpatialIndex, SpacetimeIndex, UpLo::Up,
      UpLo::Lo, UpLo::Up, ti_D, ti_j, ti_B>();

  // Rank 3: DataVector; first and second indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_ab_symmetry<
      DataVector, SpacetimeIndex, SpacetimeIndex, UpLo::Lo, UpLo::Up, ti_b,
      ti_a, ti_C>();

  // Rank 3: DataVector; first and third indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_ac_symmetry<
      DataVector, SpatialIndex, SpacetimeIndex, UpLo::Lo, UpLo::Lo, ti_i, ti_f,
      ti_j>();

  // Rank 3: DataVector; second and third indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_bc_symmetry<
      DataVector, SpacetimeIndex, SpatialIndex, UpLo::Lo, UpLo::Up, ti_d, ti_J,
      ti_I>();

  // Rank 3: DataVector; symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_3_abc_symmetry<
      DataVector, SpacetimeIndex, UpLo::Lo, ti_f, ti_d, ti_a>();

  // Rank 4: double; nonsymmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_4<
      double, Symmetry<4, 3, 2, 1>,
      index_list<SpacetimeIndex<3, UpLo::Lo, Frame::Inertial>,
                 SpacetimeIndex<3, UpLo::Up, Frame::Inertial>,
                 SpatialIndex<1, UpLo::Lo, Frame::Inertial>,
                 SpatialIndex<2, UpLo::Lo, Frame::Inertial>>,
      ti_b, ti_A, ti_k, ti_l>();

  // Rank 4: double; second and third indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_4<
      double, Symmetry<3, 2, 2, 1>,
      index_list<SpacetimeIndex<2, UpLo::Up, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Lo, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Lo, Frame::Grid>,
                 SpatialIndex<1, UpLo::Lo, Frame::Grid>>,
      ti_G, ti_d, ti_a, ti_j>();

  // Rank 4: double; first, second, and fourth indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_4<
      double, Symmetry<2, 2, 1, 2>,
      index_list<SpatialIndex<3, UpLo::Lo, Frame::Inertial>,
                 SpatialIndex<3, UpLo::Lo, Frame::Inertial>,
                 SpatialIndex<3, UpLo::Lo, Frame::Inertial>,
                 SpatialIndex<3, UpLo::Lo, Frame::Inertial>>,
      ti_j, ti_i, ti_k, ti_l>();

  // Rank 4: double; symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_4<
      double, Symmetry<1, 1, 1, 1>,
      index_list<SpacetimeIndex<3, UpLo::Up, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Up, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Up, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Up, Frame::Grid>>,
      ti_F, ti_A, ti_C, ti_D>();

  // Rank 4: DataVector; nonsymmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_4<
      DataVector, Symmetry<4, 3, 2, 1>,
      index_list<SpacetimeIndex<3, UpLo::Lo, Frame::Inertial>,
                 SpacetimeIndex<3, UpLo::Up, Frame::Inertial>,
                 SpatialIndex<1, UpLo::Lo, Frame::Inertial>,
                 SpatialIndex<2, UpLo::Lo, Frame::Inertial>>,
      ti_b, ti_A, ti_k, ti_l>();

  // Rank 4: DataVector; second and third indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_4<
      DataVector, Symmetry<3, 2, 2, 1>,
      index_list<SpacetimeIndex<2, UpLo::Up, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Lo, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Lo, Frame::Grid>,
                 SpatialIndex<1, UpLo::Lo, Frame::Grid>>,
      ti_G, ti_d, ti_a, ti_j>();

  // Rank 4: DataVector; first, second, and fourth indices symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_4<
      DataVector, Symmetry<2, 2, 1, 2>,
      index_list<SpatialIndex<3, UpLo::Lo, Frame::Inertial>,
                 SpatialIndex<3, UpLo::Lo, Frame::Inertial>,
                 SpatialIndex<3, UpLo::Lo, Frame::Inertial>,
                 SpatialIndex<3, UpLo::Lo, Frame::Inertial>>,
      ti_j, ti_i, ti_k, ti_l>();

  // Rank 4: DataVector; symmetric
  TestHelpers::TensorExpressions::test_evaluate_rank_4<
      DataVector, Symmetry<1, 1, 1, 1>,
      index_list<SpacetimeIndex<3, UpLo::Up, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Up, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Up, Frame::Grid>,
                 SpacetimeIndex<3, UpLo::Up, Frame::Grid>>,
      ti_F, ti_A, ti_C, ti_D>();
}
