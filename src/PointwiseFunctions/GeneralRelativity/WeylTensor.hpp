// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace gsl {
template <typename>
struct not_null;
}  // namespace gsl
/// \endcond

namespace gr {
/// @{
/*!
 * \ingroup GeneralRelativityGroup
 * \brief Compute the Weyl tensor from its electric and magnetic components.
 * \details We use the formula.
 *
 * We compute Weyl tensor as
 * \f[
 * C_{abcd} = \mathrm{P}
 * \Big[4 E_{ac} \left( \gamma_{bd} + n_{b} n_{d} \right) \vphantom{]}
 * \vphantom{[]} - {\epsilon_{ab}}^{e} n_{d} B_{ce} \Big]},
 * \f]
 *
 * where \f$ \mathrm{P} \f$ is a projector that imposes the symmetries of the
 * Riemann tensor --i.e. \f$R_{abcd} = R_{[ab][cd]} = R_{cdab}\f$, and
 * \f$ E_{ab} \f$ and \f$ B_{ab} \f$ are the electric and magnetic parts of the
 * Weyl tensor, respectively.
 */
template <size_t Dim, typename Frame, typename DataType>
void weyl_tensor_from_electric_magnetic(
    gsl::not_null<tnsr::abcc<DataType, Dim, Frame>*> weyl);

template <size_t Dim, typename Frame, typename DataType>
tnsr::abcc<DataType, Dim, Frame> weyl_tensor_from_electric_magnetic();
/// @}

/// @{
/*!
 * \ingroup GeneralRelativityGroup
 * \brief Compute the square of the Weyl tensor.
 * \details We use the formula.
 *
 * We compute the squared Weyl tensor \f$ C_{abcd} C^{abcd} \f^$.
 */
template <typename DataType, typename Frame>
void weyl_square_scalar_compute(gsl::not_null<Scalar<DataType>*> weyl_square);

template <typename DataType, typename Frame>
Scalar<DataType> weyl_square_scalar_compute();
/// @}

namespace Tags {

/// Compute item for the Weyl tensor computed from the electric and magnetic
/// parts of the Weyl tensor in vacuum.
///
/// Can be retrieved using gr::Tags::Weyl
template <size_t Dim, typename Frame, typename DataType>
struct WeylCompute : Weyl<Dim, Frame, DataType>, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::WeylMagneticCompute<Frame, DataType>,
                 gr::Tags::WeylElectricCompute<Frame, DataType>,
                 gr::Tags::SpatialMetric<SpatialDim, Frame, DataType>>;

  using return_type = tnsr::abcc<DataType, Dim, Frame>;

  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<tnsr::abcc<DataType, Dim, Frame>*>
      /* more input arguments */)>(
          &weyl_tensor_from_electric_magnetic<Dim, Frame, DataType>);

  using base = Weyl<Dim, Frame, DataType>;
};

/// Compute item for the Weyl squared scalar from the Weyl tensor.
///
/// Can be retrieved using gr::Tags::WeylSquareScalar
template <typename DataType, typename Frame>
struct WeylSquareScalarCompute : WeylSquareScalar<DataType>, db::ComputeTag {
  using argument_tags =
      tmpl::list<gr::Tags::WeylCompute<4_st, Frame, DataType>,
                 gr::Tags::InverseSpatialMetric<3_st, Frame, DataType>>;

  using return_type = Scalar<DataType>;

  static constexpr auto function = static_cast<void (*)(
      gsl::not_null<Scalar<DataType>*> /* more input arguments */)>(
      &weyl_square_scalar_compute<DataType, Frame>);

  using base = WeylSquareScalar<DataType>;
};

}  // namespace gr
