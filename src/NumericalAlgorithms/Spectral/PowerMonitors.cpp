// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "NumericalAlgorithms/Spectral/PowerMonitors.hpp"

#include <cmath>

#include "Parallel/Printf.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace PowerMonitors {

void compute_power_monitor( const gsl::not_null<double*> result) {
    Parallel::printf("Inside compute monitor function");
    *result = 1.0;
}

} // namespace PowerMonitors
