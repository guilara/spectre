// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>
#include <string>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/FixedScalarTensor/FixedDecoupledScalar/Solutions.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrSchild.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma:  no_include <pup.h>

/// \cond
namespace PUP {
class er;  // IWYU pragma: keep
}  // namespace PUP
/// \endcond

namespace fe::DecoupledScalar::Solutions {
/*!
 * \brief Analytic initial data for a pure spherical harmonic in three
 * dimensions.
 *
 * \details The initial data is taken from \cite Scheel2003vs , Eqs. 4.1--4.3,
 * and sets the evolved variables of the scalar wave as follows:
 *
 * \f{align}
 * \Psi &= 0 \\
 * \Phi_i &= 0 \\
 * \Pi &= \Pi_0(r, \theta, \phi) =  A e^{- (r - r_0)^2 / w^2} Y_{lm}(\theta,
 * \phi), \f}
 *
 * where \f$r_0\f$ is the radius of the profile and \f$w\f$ is its width. This
 * describes a pure spherical harmonic mode \f$Y_{lm}(\theta, \phi)\f$ truncated
 * by a circular Gaussian window function.
 *
 * When evolved, the scalar field \f$\Phi\f$ will briefly build up around the
 * radius \f$r_0\f$ and then disperse. This can be used to study the ringdown
 * behavior and late-time tails in different background spacetimes.
 */
class KerrSphericalHarmonic
    : /* public evolution::initial_data::InitialData, */
      /* Why does it work with `virtual`? */
      public virtual evolution::initial_data::InitialData,
      public AnalyticSolution,
      public MarkAsAnalyticSolution {
 public:
  /// The mass of the black hole.
  struct Mass {
    using type = double;
    static constexpr Options::String help = {"Mass of the black hole."};
    static type lower_bound() { return 0.0; }
  };
  /// The spin of the black hole
  struct Spin {
    using type = std::array<double, volume_dim>;
    static constexpr Options::String help = {
        "The [x,y,z] dimensionless spin of the black hole"};
  };
  /// The amplitude of the scalar field
  struct Amplitude {
    using type = double;
    static constexpr Options::String help = {
        "Amplitude of the constant scalar field"};
  };
  /// The amplitude of the scalar driver
  struct AmplitudeOfDriver {
    using type = double;
    static constexpr Options::String help = {
        "Amplitude of the scalar driver field"};
  };
  struct Radius {
    using type = double;
    static constexpr Options::String help = {
        "The radius of the spherical harmonic profile"};
    static type lower_bound() { return 0.0; }
  };

  struct Width {
    using type = double;
    static constexpr Options::String help = {
        "The width of the spherical harmonic profile. The width must be "
        "greater than 0."};
    static type lower_bound() { return 0.0; }
  };

  struct Mode {
    using type = std::pair<size_t, int>;
    static constexpr Options::String help = {
        "The l-mode and m-mode of the spherical harmonic Ylm. The absolute "
        "value of the m_mode must be less than or equal to the "
        "l-mode."};
  };

  using options =
      tmpl::list<Mass, Spin, Amplitude, AmplitudeOfDriver, Radius, Width, Mode>;
  static constexpr Options::String help = {
      "Initial data for a pure spherical harmonic mode truncated by a circular "
      "Gaussian window funtion. The expression is taken from Scheel(2003), "
      "equations 4.1-4.3."};

  KerrSphericalHarmonic() = default;
  KerrSphericalHarmonic(const KerrSphericalHarmonic& /*rhs*/) = default;
  KerrSphericalHarmonic& operator=(const KerrSphericalHarmonic& /*rhs*/) =
      default;
  KerrSphericalHarmonic(KerrSphericalHarmonic&& /*rhs*/) = default;
  KerrSphericalHarmonic& operator=(KerrSphericalHarmonic&& /*rhs*/) = default;
  ~KerrSphericalHarmonic() override = default;

  KerrSphericalHarmonic(double mass,
                        const std::array<double, 3>& dimensionless_spin,
                        double amplitude, double amplitude_of_driver,
                        double radius, double width,
                        std::pair<size_t, int> mode);

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit KerrSphericalHarmonic(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(KerrSphericalHarmonic);
  /// \endcond

  // The extra tags below can also be added as common to all scalar tensor
  // solutions
  template <typename DataType, typename Frame = Frame::Inertial>
  using tags = tmpl::flatten<tmpl::list<
      typename AnalyticSolution::template tags<DataType>,
      gr::Tags::DerivDetSpatialMetric<DataType, 3, Frame>,
      gr::Tags::TraceExtrinsicCurvature<DataType>,
      gr::Tags::SpatialChristoffelFirstKind<DataType, 3, Frame>,
      gr::Tags::SpatialChristoffelSecondKind<DataType, 3, Frame>,
      gr::Tags::TraceSpatialChristoffelSecondKind<DataType, 3, Frame>>>;

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

  /// @{
  /// Retrieve scalar driver variable at `(x, t)`
  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<fe::ScalarDriver::Tags::Psi> /*meta*/) const
      -> tuples::TaggedTuple<fe::ScalarDriver::Tags::Psi>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<fe::ScalarDriver::Tags::Phi<3_st>> /*meta*/) const
      -> tuples::TaggedTuple<fe::ScalarDriver::Tags::Phi<3_st>>;

  template <typename DataType>
  auto variables(const tnsr::I<DataType, 3>& x, double t,
                 tmpl::list<fe::ScalarDriver::Tags::Pi> /*meta*/) const
      -> tuples::TaggedTuple<fe::ScalarDriver::Tags::Pi>;
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
  friend bool operator==(const KerrSphericalHarmonic& lhs,
                         const KerrSphericalHarmonic& rhs);

  double mass_ = std::numeric_limits<double>::signaling_NaN();
  std::array<double, volume_dim> dimensionless_spin_ =
      make_array<volume_dim>(std::numeric_limits<double>::signaling_NaN());
  double amplitude_ = std::numeric_limits<double>::signaling_NaN();
  double amplitude_of_driver_ = std::numeric_limits<double>::signaling_NaN();
  double radius_{std::numeric_limits<double>::signaling_NaN()};
  double width_sq_{std::numeric_limits<double>::signaling_NaN()};
  std::pair<size_t, int> mode_{std::numeric_limits<size_t>::signaling_NaN(),
                               std::numeric_limits<int>::signaling_NaN()};
  // Add BH parameters
  gr::Solutions::KerrSchild background_spacetime_{};
};

bool operator!=(const KerrSphericalHarmonic& lhs,
                const KerrSphericalHarmonic& rhs);

}  // namespace fe::DecoupledScalar::Solutions
