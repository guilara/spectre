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
 * \brief Returns array of power monitors in each dimension.
 *
 * The are computed following Ref. \cite Szilagyi2014fna, e.g. in the
 * x dimension, we compute
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
std::array<DataVector, Dim> power_monitor_array(const DataVector&,
                                                const Mesh<Dim>&);
/// @}

}  // namespace PowerMonitors
