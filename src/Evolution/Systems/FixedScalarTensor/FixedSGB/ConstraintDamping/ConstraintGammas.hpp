// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/ConstraintDamping/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/FixedSGB/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

// IWYU pragma: no_forward_declare Tags::deriv

/// \cond
namespace Tags {
struct Time;
}  // namespace Tags
namespace domain {
namespace Tags {
template <size_t Dim, typename Frame>
struct Coordinates;
}  // namespace Tags
}  // namespace domain
class DataVector;
template <typename X, typename Symm, typename IndexList>
class Tensor;
/// \endcond

namespace fe::sgb::ConstraintDamping::Tags {

/*!
 * \brief Computes the constraint damping parameter \f$\sigma\f$ from the
 * coordinates and a DampingFunction.
 *
 * \details Can be retrieved using
 * `fe::ScalarTensorDriver::Tags::SigmaParameter`.
 */
template <size_t SpatialDim, typename Frame>
struct SigmaParameterCompute : ::fe::ScalarTensorDriver::Tags::SigmaParameter,
                               db::ComputeTag {
  using argument_tags =
      tmpl::list<DampingFunctionSigmaParameter<SpatialDim, Frame>,
                 domain::Tags::Coordinates<SpatialDim, Frame>, ::Tags::Time,
                 ::domain::Tags::FunctionsOfTime>;
  using return_type = Scalar<DataVector>;

  static constexpr void function(
      const gsl::not_null<Scalar<DataVector>*> gamma,
      const ::fe::sgb::ConstraintDamping::DampingFunction<SpatialDim, Frame>&
          damping_function,
      const tnsr::I<DataVector, SpatialDim, Frame>& coords, const double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) {
    damping_function(gamma, coords, time, functions_of_time);
  }

  using base = ::fe::ScalarTensorDriver::Tags::SigmaParameter;
};

/*!
 * \brief Computes the constraint damping parameter \f$\tau\f$ from the
 * coordinates and a DampingFunction.
 *
 * \details Can be retrieved using
 * `fe::ScalarTensorDriver::Tags::TauParameter`.
 */
template <size_t SpatialDim, typename Frame>
struct TauParameterCompute : ::fe::ScalarTensorDriver::Tags::TauParameter,
                             db::ComputeTag {
  using argument_tags =
      tmpl::list<DampingFunctionTauParameter<SpatialDim, Frame>,
                 domain::Tags::Coordinates<SpatialDim, Frame>, ::Tags::Time,
                 ::domain::Tags::FunctionsOfTime>;
  using return_type = Scalar<DataVector>;

  static constexpr void function(
      const gsl::not_null<Scalar<DataVector>*> gamma,
      const ::fe::sgb::ConstraintDamping::DampingFunction<SpatialDim, Frame>&
          damping_function,
      const tnsr::I<DataVector, SpatialDim, Frame>& coords, const double time,
      const std::unordered_map<
          std::string,
          std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
          functions_of_time) {
    damping_function(gamma, coords, time, functions_of_time);
  }

  using base = ::fe::ScalarTensorDriver::Tags::TauParameter;
};

}  // namespace fe::sgb::ConstraintDamping::Tags
