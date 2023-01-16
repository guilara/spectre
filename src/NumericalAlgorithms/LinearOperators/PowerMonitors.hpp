// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
/// \endcond

namespace PowerMonitors {

/// @{
/*!
 * \ingroup NumericalAlgorithmsGroup
 * \brief Computes the array of power monitors from DataBox variables.
 *
 */
template <size_t Dim>
void compute_power_monitor(gsl::not_null<Scalar<DataVector>*> result,
                           const Scalar<DataVector>&,
                           const tnsr::i<DataVector, Dim, Frame::Inertial>&,
                           const Mesh<Dim>&);
/// @}

/// @{
/*!
 * \ingroup NumericalAlgorithmsGroup
 * \brief Returns array of power monitors in each dimension.
 *
 * The are computed following Ref. \cite{Szilagyi2014fna}, e.g. in the
 * x dimension, we compute
 *
 * \f{align*}{
 *  \log10 P_{k_0}[\psi] = \log10 \sqrt{ \frac{1}{N_1 N_2}
 *   \sum_{k_1,k_2} \left| C_{k_0,k_1,k_2} \right|^2} ,
 * \f}
 *
 * where \f$ C_{k_0,k_1,k_2}\f$ are the modal coefficients
 * of variable \f$ \psi \f$.
 *
 */
template <size_t Dim>
std::array<DataVector, Dim> power_monitor_array(const DataVector&,
                                                const Mesh<Dim>&);
/// @}

namespace Tags {
template <size_t Dim>
struct PowerMonitor : db::SimpleTag {
  using type = Scalar<DataVector>;
};

/// @{
/*!
 * \ingroup NumericalAlgorithmsGroup
 * \brief Compute tag for power monitors.
 *
 */
template <size_t Dim>
struct PowerMonitorCompute : PowerMonitor<Dim>, db::ComputeTag {
  using argument_tags =
      tmpl::list<ScalarWave::Tags::Pi, ScalarWave::Tags::Phi<Dim>,
                 domain::Tags::Mesh<Dim>>;
  using base = PowerMonitor<Dim>;
  using return_type = Scalar<DataVector>;

  static constexpr void (*function)(
      const gsl::not_null<return_type*> result, const Scalar<DataVector>&,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&,
      const Mesh<Dim>&) = &compute_power_monitor<Dim>;
};
/// @}

}  // namespace Tags
}  // namespace PowerMonitors
