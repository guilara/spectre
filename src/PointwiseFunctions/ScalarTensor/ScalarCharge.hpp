// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <utility>

#include "ApparentHorizons/StrahlkorperGr.hpp"
#include "ApparentHorizons/TagsDeclarations.hpp"  // IWYU pragma: keep
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Strahlkorper.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/TagsDeclarations.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/TagsTypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/TagsDeclarations.hpp"  // IWYU pragma: keep
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace ScalarTensor {

/// @{
/*!
 * \brief Computes the dot product between the scalar gradient and unit normal
 * the surface.
 *
 * \details Compute the integrand of:
 * \f{align*}
 * q = - \dfrac{1}{4 \pi} \oint dA \Phi_i n^{i},
 * \f}
 * where \f$ n^{i} \f$ is the unit (outward?) normal of the surface.
 *
 * This makes sense if \f$ \psi \sim \phi_\infty + q / r + \cdots \f$, for which
 * \f$ \partial_r \psi \sim  - q / r^2 + \cdots \f$, and thus
 * \f$ q \sim - r^2 \partial_r \psi  \f$.
 *
 */
template <typename Frame>
void scalar_charge_integrand(
    const gsl::not_null<Scalar<DataVector>*> result,
    const tnsr::i<DataVector, 3, Frame>& phi,
    const tnsr::I<DataVector, 3, Frame>& unit_normal_vector,
    const Scalar<DataVector>& area_element);
/// @}

} // namespace ScalarTensor

namespace ScalarTensor::StrahlkorperScalar::Tags {

/// Scalar charge tag.
struct ScalarChargeIntegrand : db::SimpleTag {
  using type = Scalar<DataVector>;
};

/// Calculates the scalar charge.
template <typename Frame>
struct ScalarChargeIntegrandCompute : ScalarChargeIntegrand, db::ComputeTag {
  using base = ScalarChargeIntegrand;
  static constexpr auto function = static_cast<void (*)(
      gsl::not_null<Scalar<DataVector>*>,
      const tnsr::i<DataVector, 3_st, Frame>&,
      const tnsr::I<DataVector, 3_st, Frame>&, const Scalar<DataVector>&)>(
      &ScalarTensor::scalar_charge_integrand<Frame>);
  using argument_tags = tmpl::list<CurvedScalarWave::Tags::Phi<3_st>,
                                   StrahlkorperTags::UnitNormalVector<Frame>,
                                   StrahlkorperGr::Tags::AreaElement<Frame>>;
  using return_type = Scalar<DataVector>;
};

}  // namespace ScalarTensor::StrahlkorperTags
