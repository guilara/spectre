// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/Spectral/PowerMonitors.hpp"

#include <cmath>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace PowerMonitors {

void compute_power_monitor(const gsl::not_null<Scalar<DataVector>*> result,
                           const Scalar<DataVector>& pi) {
    Parallel::printf("Inside compute monitor function \n");
    Parallel::printf("get_size(pi) = %u \n", get_size(get(pi)));
    destructive_resize_components(result, get_size(get(pi)));
    get(*result) += square(get(pi));
}

} // namespace PowerMonitors
