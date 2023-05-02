// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>
#include <string>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/GaugeWave.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/Solutions.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma:  no_include <pup.h>

/// \cond
namespace PUP {
class er;  // IWYU pragma: keep
}  // namespace PUP
/// \endcond

namespace ScalarTensor::Solutions {
/*!
 * \brief Set the scalar variables to zero on Minkowski in coordinates with a
 * Gauge Wave.
 */
class GaugeWaveConstantScalar
    : /* public evolution::initial_data::InitialData, */
      /* Why does it work with `virtual`? */
      public virtual evolution::initial_data::InitialData,
      public AnalyticSolution,
      public MarkAsAnalyticSolution {
 public:
  /// The amplitude of the gauge wave
  struct AmplitudeGaugeWave {
    using type = double;
    static constexpr Options::String help = {
        "The amplitude of the gauge wave."};
    static type upper_bound() { return 1.; }
    static type lower_bound() { return -1.; }
  };
  /// The wavelength of the gauge wave
  struct WavelengthGaugeWave {
    using type = double;
    static constexpr Options::String help = {
        "The wavelength of the gauge wave."};
    static type lower_bound() { return 0.; }
  };
  /// The amplitude of the scalar field
  struct AmplitudeScalar {
    using type = double;
    static constexpr Options::String help = {
        "The amplitude of the scalar field."};
  };

  using options =
      tmpl::list<AmplitudeGaugeWave, WavelengthGaugeWave, AmplitudeScalar>;
  static constexpr Options::String help = {
      "Zero scalar field in Minkowski in Gauge Wave coordinates."};

  GaugeWaveConstantScalar() = default;
  GaugeWaveConstantScalar(const GaugeWaveConstantScalar& /*rhs*/) = default;
  GaugeWaveConstantScalar& operator=(const GaugeWaveConstantScalar& /*rhs*/) =
      default;
  GaugeWaveConstantScalar(GaugeWaveConstantScalar&& /*rhs*/) = default;
  GaugeWaveConstantScalar& operator=(GaugeWaveConstantScalar&& /*rhs*/) =
      default;
  ~GaugeWaveConstantScalar() override = default;

  GaugeWaveConstantScalar(double amplitude_gauge_wave,
                          double wavelength_gauge_wave,
                          double amplitude_scalar);

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit GaugeWaveConstantScalar(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(GaugeWaveConstantScalar);
  /// \endcond

  // The extra tags below can also be added as common to all scalar tensor
  // solutions
  template <typename DataType, typename Frame = Frame::Inertial>
  using tags = tmpl::flatten<tmpl::list<
      typename AnalyticSolution::template tags<DataType>
      // The computations for the following tags need to be added
    //   ,
    //   gr::Tags::DerivDetSpatialMetric<3_st, Frame, DataType>,
    //   gr::Tags::TraceExtrinsicCurvature<DataType>,
    //   gr::Tags::SpatialChristoffelFirstKind<3_st, Frame, DataType>,
    //   gr::Tags::SpatialChristoffelSecondKind<3_st, Frame, DataType>,
    //   gr::Tags::TraceSpatialChristoffelSecondKind<3_st, Frame, DataType>
      >>;

  /// @{
  /// Retrieve scalar variable at `(x, t)`
  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<CurvedScalarWave::Tags::Psi> /*meta*/) const
      -> tuples::TaggedTuple<CurvedScalarWave::Tags::Psi>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<CurvedScalarWave::Tags::Phi<3_st>> /*meta*/) const
      -> tuples::TaggedTuple<CurvedScalarWave::Tags::Phi<3_st>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<CurvedScalarWave::Tags::Pi> /*meta*/) const
      -> tuples::TaggedTuple<CurvedScalarWave::Tags::Pi>;
  /// @}

  /// Retrieve a collection of scalar variables at `(x, t)`
  template <typename DataType, typename... Tags>
  tuples::TaggedTuple<Tags...> variables(const tnsr::I<DataType, 3>& x,
                                         double t,
                                         tmpl::list<Tags...> /*meta*/) const {
    static_assert(sizeof...(Tags) > 1,
                  "The generic template will recurse infinitely if only one "
                  "tag is being retrieved.");
    return {get<Tags>(variables(x, t, tmpl::list<Tags>{}))...};
  }

  /// Retrieve the metric variables
  template <typename DataType, typename Tag>
  tuples::TaggedTuple<Tag> variables(const tnsr::I<DataType, 3>& x, double t,
                                     tmpl::list<Tag> /*meta*/) const {
    return background_spacetime_.variables(x, t, tmpl::list<Tag>{});
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) override;

 protected:
  friend bool operator==(const GaugeWaveConstantScalar& lhs,
                         const GaugeWaveConstantScalar& rhs);

  double amplitude_gauge_wave_ = std::numeric_limits<double>::signaling_NaN();
  double wavelength_gauge_wave_ = std::numeric_limits<double>::signaling_NaN();
  double amplitude_scalar_ = std::numeric_limits<double>::signaling_NaN();
  gr::Solutions::GaugeWave<3> background_spacetime_{};
};

bool operator!=(const GaugeWaveConstantScalar& lhs,
                const GaugeWaveConstantScalar& rhs);

}  // namespace ScalarTensor::Solutions
