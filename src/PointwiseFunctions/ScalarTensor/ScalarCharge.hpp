// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

namespace ScalarTensor {

template <typename Frame>
void scalar_charge_integrand(
    const gsl::not_null<Scalar<DataVector>*> result,
    const tnsr::i<DataVector, 3, Frame>& phi,
    const tnsr::I<DataVector, 3, Frame>& unit_normal_vector,
    const Scalar<DataVector>& area_element);

} // namespace ScalarTensor

namespace ScalarTensor::StrahlkorperScalar::Tags {

/// Scalar charge tag.
struct ScalarCharge : db::SimpleTag {
  using type = Scalar<DataVector>;
};

/// Calculates the scalar charge.
template <typename Frame>
struct ScalarChargeCompute : ScalarCharge, db::ComputeTag {
  using base = ScalarCharge;
  static constexpr auto function = static_cast<void (*)(
      gsl::not_null<Scalar<DataVector>*>, const tnsr::i<DataVector, 3, Frame>&,
      const tnsr::I<DataVector, 3, Frame>&, const Scalar<DataVector>&)>(
      &StrahlkorperGr::scalar_charge<Frame>);
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Phi<3_st>,
                 StrahlkorperTags::UnitNormalVector<Frame>, AreaElement<Frame>>;
  using return_type = Scalar<DataVector>;
};

}  // namespace ScalarTensor::StrahlkorperTags
