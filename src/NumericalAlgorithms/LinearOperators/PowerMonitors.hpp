// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
/// \endcond

/*!
 * \brief Items for assessing truncation error in spectral methods.
 */
namespace PowerMonitors {

namespace detail {

double relative_truncation_error_impl(const DataVector& input_power_monitors,
                                      const size_t upperBound);
} // namespace detail

/// @{
/*!
 * \ingroup SpectralGroup
 * \brief Returns array of power monitors in each spatial dimension.
 *
 * Computed following Sec. 5.1 of Ref. \cite Szilagyi2014fna.
 * For example, in the x dimension (indexed by \f$ k_0 \f$), we compute
 *
 * \f{align*}{
 *  P_{k_0}[\psi] = \sqrt{ \frac{1}{N_1 N_2}
 *   \sum_{k_1,k_2} \left| C_{k_0,k_1,k_2} \right|^2} ,
 * \f}
 *
 * where \f$ C_{k_0,k_1,k_2}\f$ are the modal coefficients
 * of variable \f$ \psi \f$.
 *
 */
template <size_t Dim>
void power_monitors(gsl::not_null<std::array<DataVector, Dim>*> result,
                    const DataVector& input_data_vector, const Mesh<Dim>& mesh);

template <size_t Dim>
std::array<DataVector, Dim> power_monitors(const DataVector& input_data_vector,
                                           const Mesh<Dim>& mesh);
/// @}

/*!
 * \ingroup SpectralGroup
 * \brief Returns the maximum value of a variable in the element.
 *
 * Compute the maximum of variable \f$ u \f$ in the element as
 *
 * \f{align*}{
 *  u_\mathrm{max} = \max_{ijk} \left\lvert u_{ijk} \right\rvert ,
 * \f}
 *
 * where \f$ u_{ijk}\f$ are the nodal coefficients
 * of variable \f$ u \f$.
 *
 */
double maximum_of_variable(const DataVector& input_data_vector);

/// @{
/*!
 * \ingroup SpectralGroup
 * \brief Relative truncation error.
 *
 * Truncation error according to Eqs. (57) and (58) of Ref. \cite
 * Szilagyi2014fna, i.e.,
 *
 * \f{align*}{
 *  \mathcal{T}\left[P_k\right] = \log_{10} \max \left(P_0, P_1\right)
 *   - \dfrac{\sum_{j} \log_{10} \left(P_j\right) w_j }{\sum_{j} w_j} ,
 * \f}
 *
 * with weights
 *
 * \f{align*}{
 *  w_j = \exp\left[ - (j - N_k + \dfrac{1}{2})^2 \right] .
 * \f}
 *
 * where \f$ N_k \f$ is the number of power monitors. Here the second term
 * is a weighted average with larger weights toward the highest modes.
 *
 */
template <size_t Dim>
void relative_truncation_error(gsl::not_null<std::array<double, Dim>*> result,
                      const DataVector& input_data_vector,
                      const Mesh<Dim>& mesh);

template <size_t Dim>
std::array<double, Dim> relative_truncation_error(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh);
/// @}

/// @{
/*!
 * \ingroup SpectralGroup
 * \brief Returns an estimate of the absolute truncation error.
 *
 * The estimate of the numerical error is given by
 *
 * \f{align*}{
 *  \mathcal{E}\left[P_k\right] = \dfrac{u_\mathrm{max} \times 10^{
 *   - \mathcal{T}[P_k]}}{
 *  \mathrm{atol} + \mathrm{rtol} \times
 *  u_\mathrm{max}} ,
 * \f}
 *
 * where \f$ \mathcal{T}[P_k] \f$ is the relative error estimate computed from
 * the power monitors \f$ P_k \f$. We set \f$ \mathrm{rtol} = 0 \f$ and
 * \f$ \mathrm{rtol} = 1 \f$.
 *
 */
template <size_t Dim>
void error_estimate(gsl::not_null<std::array<double, Dim>*> result,
                               const DataVector& input_data_vector,
                               const Mesh<Dim>& mesh);

template <size_t Dim>
std::array<double, Dim> error_estimate(const DataVector& input_data_vector,
                                       const Mesh<Dim>& mesh);
/// @}

/// @{
/*!
 * \ingroup SpectralGroup
 * \brief Computes the minimum of two absolute truncation error estimates.
 *
 * Calculated as smallest relative truncation error for all power monitors and
 * with one less
 *
 * \f{align*}{
 *  \epsilon\left[ P_k \right] = \min \left(\mathcal{T}_{N_k}\left[P_k\right],
 * \mathcal{T}_{N_k - 1}\left[P_k\right]\right) ,
 * \f}
 *
 * where \f$ N_k \f$ is the number of power monitors.
 */
template <size_t Dim>
void absolute_truncation_error_estimate(
    gsl::not_null<std::array<double, Dim>*> result,
    const DataVector& input_data_vector, const Mesh<Dim>& mesh);

template <size_t Dim>
std::array<double, Dim> absolute_truncation_error_estimate(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh);
/// @}

}  // namespace PowerMonitors
