// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <pup.h>
#include <utility>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Options/Options.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace CurvedScalarWave::AnalyticData {

/*!
 * \brief Analytic initial data for scalar that is zero except for random noise.
 *
 * \details The initial data is for now a constant.
 */

class NoiseFloorScalar : public MarkAsAnalyticData {
 public:
  struct Amplitude {
    using type = double;
    static constexpr Options::String help = {
        "The amplitude of the random noise fluctuation"};
    static type lower_bound() { return 0.0; }
  };

  using options = tmpl::list<Amplitude>;

  static constexpr Options::String help = {
      "Initial data for a scalar with random noise values."};

  NoiseFloorScalar() = default;

  NoiseFloorScalar(double amplitude, const Options::Context& context = {});

  static constexpr size_t volume_dim = 3;
  using tags =
      tmpl::list<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
                 CurvedScalarWave::Tags::Phi<3>>;

  /// Retrieve the evolution variables at spatial coordinates `x`
  tuples::TaggedTuple<CurvedScalarWave::Tags::Psi, CurvedScalarWave::Tags::Pi,
                      CurvedScalarWave::Tags::Phi<3>>
  variables(const tnsr::I<DataVector, 3>& x, tags /*meta*/) const;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/);

 private:
  double amplitude_{std::numeric_limits<double>::signaling_NaN()};

  friend bool operator==(const NoiseFloorScalar& lhs,
                         const NoiseFloorScalar& rhs);

  friend bool operator!=(const NoiseFloorScalar& lhs,
                         const NoiseFloorScalar& rhs);
};

}  // namespace CurvedScalarWave::AnalyticData
