// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarTensor/Actions/NumericInitialData.hpp"

#include <boost/functional/hash.hpp>
#include <string>
#include <utility>
#include <variant>

namespace ScalarTensor {

NumericInitialData::NumericInitialData(
    std::string file_glob, std::string subfile_name,
    std::variant<double, importers::ObservationSelector> observation_value,
    const bool enable_interpolation,
    typename GhNumericId::Variables::type gh_selected_variables,
    // typename HydroNumericId::Variables::type hydro_selected_variables,
    typename ScalarVariables scalar_selected_variables,
    const double density_cutoff)
    : gh_numeric_id_(file_glob, subfile_name, observation_value,
                     enable_interpolation, std::move(gh_selected_variables)),
      scalar_selected_variables_(std::move(scalar_selected_variables)) {}

NumericInitialData::NumericInitialData(CkMigrateMessage* msg)
    : InitialData(msg) {}

PUP::able::PUP_ID NumericInitialData::my_PUP_ID = 0;

size_t NumericInitialData::volume_data_id() const {
  size_t hash = 0;
  boost::hash_combine(hash, gh_numeric_id_.volume_data_id());
  boost::hash_combine(hash, "ScalarID");
  return hash;
}

void NumericInitialData::pup(PUP::er& p) {
  p | gh_numeric_id_;
  p | scalar_selected_variables_;
//   p | hydro_numeric_id_;
}

bool operator==(const NumericInitialData& lhs, const NumericInitialData& rhs) {
  return lhs.gh_numeric_id_ == rhs.gh_numeric_id_ and
        //  lhs.hydro_numeric_id_ == rhs.hydro_numeric_id_;
        lhs.scalar_selected_variables_ == rhs.scalar_selected_variables_;
}

}  // namespace ScalarTensor
