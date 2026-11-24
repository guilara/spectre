// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/Tensor/TypeAliases.hpp"

namespace sr {
/*!
 * \ingroup SpecialRelativityGroup
 * \brief Computes the matrix that undoes the Lorentz contraction of an object
 * moving with constant velocity \f$v^i\f$.
 *
 * \details Consider an object at rest in coordinates \f$(T, X^i)\f$ and
 * moving with velocity \f$v^i\f$ (with \f$c=1\f$) in coordinates
 * \f$(t, x^i)\f$ related by a Lorentz boost. Let \f$r^i = x^i - v^i t\f$ be
 * the position relative to the moving center on a slice of constant \f$t\f$.
 * Then the rest-frame position is \f$X^i = \Lambda^i{}_j r^j\f$ with
 *
 * \f{align}{
 * \Lambda^i{}_j = \delta^i_j + (\gamma - 1) \hat{n}^i \hat{n}_j
 *   = \delta^i_j + \frac{\gamma^2}{1 + \gamma} v^i v_j,
 * \f}
 *
 * where \f$\hat{n}^i = v^i / v\f$, \f$v = \sqrt{\delta_{ij} v^i v^j}\f$ and
 * \f$\gamma = 1/\sqrt{1-v^2}\f$. This is the spatial block
 * \f$\Lambda^i{}_{\bar{j}}\f$ of `sr::lorentz_boost_matrix`, which is used to
 * compute it. The second form avoids the division by \f$v^2\f$, which loses
 * precision for small velocities.
 *
 * The contraction matrix stretches lengths along \f$\hat{n}^i\f$ by
 * \f$\gamma\f$ and leaves lengths orthogonal to \f$\hat{n}^i\f$ unchanged. Its
 * inverse, which applies the Lorentz contraction, is computed by
 * `sr::inverse_lorentz_contraction_matrix`.
 */
template <size_t SpatialDim>
tnsr::Ij<double, SpatialDim, Frame::NoFrame> lorentz_contraction_matrix(
    const std::array<double, SpatialDim>& velocity);

/*!
 * \ingroup SpecialRelativityGroup
 * \brief Computes the inverse of `sr::lorentz_contraction_matrix`, which
 * applies the Lorentz contraction of an object moving with constant velocity
 * \f$v^i\f$.
 *
 * \details The inverse contraction matrix maps the rest-frame position
 * \f$X^i\f$ to the position \f$r^i = (\Lambda^{-1})^i{}_j X^j\f$ relative to
 * the moving center (see `sr::lorentz_contraction_matrix`):
 *
 * \f{align}{
 * (\Lambda^{-1})^i{}_j = \delta^i_j - \left(1 - \gamma^{-1}\right)
 *   \hat{n}^i \hat{n}_j
 *   = \delta^i_j - \frac{(\gamma v^i)(\gamma v_j)}{\gamma (1 + \gamma)}.
 * \f}
 *
 * The second form is evaluated with the components
 * \f$\Lambda^t{}_{\bar{t}} = \gamma\f$ and
 * \f$\Lambda^i{}_{\bar{t}} = \gamma v^i\f$ of `sr::lorentz_boost_matrix`. It
 * avoids both the division by \f$v^2\f$ and the cancellation in
 * \f$1 - \gamma^{-1}\f$ for small velocities.
 */
template <size_t SpatialDim>
tnsr::Ij<double, SpatialDim, Frame::NoFrame> inverse_lorentz_contraction_matrix(
    const std::array<double, SpatialDim>& velocity);
}  // namespace sr
