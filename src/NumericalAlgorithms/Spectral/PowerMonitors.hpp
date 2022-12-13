// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
// #include "Domain/Tags.hpp"
// #include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

// Implement a function that computes the power monitors
// Take as input a DataVector and a Mesh
// Return std::array<DataVector, Dim> withe the power monitors
// in each dimension

// Overload with input Variables

// Expose the function to Python for use in post-processing scripts

// Add a compute item.
// Add it to the DataBox to compute the power monitors during a simulation
// Output the power monitors to H5 files as 1D volume data
// Add Python script to read the H5 files and plot the monitors

/// \cond
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
/// \endcond

namespace PowerMonitors {

// template <size_t Dim>
// size_t compute_power_monitor (Mesh<Dim>&);

void compute_power_monitor(gsl::not_null<double*> result);

namespace Tags {

struct PowerMonitor : db::SimpleTag {
//   using type = size_t;
    using type = double;
};

// template <size_t Dim>
struct PowerMonitorCompute : PowerMonitor, db::ComputeTag {
    using argument_tags = tmpl::list<>;
    using base = PowerMonitor;
    using return_type = base::type;

    static constexpr void (*function) (
        const gsl::not_null<return_type*> result) = &compute_power_monitor ;
};

} // namespace Tags
}  // namespace PowerMonitors
