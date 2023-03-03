// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"
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

/// @{
/*!
 * \ingroup SpectralGroup
 * \brief Compute the negative log10 of the relative truncation error.
 *
 * Truncation error according to Eqs. (57) and (58) of Ref.
 * \cite Szilagyi2014fna, i.e.,
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
 * where \f$ N_k \f$ is the number of modes or gridpoints in dimension k. Here
 * the second term is a weighted average with larger weights toward the highest
 * modes. The negative of this number should correspond to the number of digits
 * resolved by the spectral expansion.
 *
 */
double relative_truncation_error(const DataVector& input_power_monitors,
                                 const size_t num_modes);

template <size_t Dim>
std::array<double, Dim> relative_truncation_error(
    const DataVector& input_data_vector, const Mesh<Dim>& mesh);

template <size_t Dim>
std::array<double, Dim> relative_truncation_error(
    const std::array<DataVector, Dim>& input_power_monitors);
/// @}

/// @{
/*!
 * \ingroup SpectralGroup
 * \brief Returns an estimate of the absolute truncation error.
 *
 * The estimate of the numerical error is given by
 *
 * \f{align*}{
 *  \mathcal{E}\left[P_k\right] = u_\mathrm{max} \times 10^{- \mathcal{T}[P_k]},
 * \f}
 *
 * where \f$ \mathcal{T}[P_k] \f$ is the relative error estimate computed from
 * the power monitors \f$ P_k \f$.
 */
template <size_t Dim>
std::array<double, Dim> truncation_error(const DataVector& input_data_vector,
                                       const Mesh<Dim>& mesh);

template <size_t Dim>
std::array<double, Dim> truncation_error(
    const std::array<double, Dim>& relative_truncation_error,
    const DataVector& input_data_vector);
/// @}

template <size_t Dim>
double truncation_error_max(const Scalar<DataVector>& input_scalar,
                            const Mesh<Dim>& mesh);

}  // namespace PowerMonitors
