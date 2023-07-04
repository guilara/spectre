// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <limits>
#include <string>

#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticSolutions/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Minkowski.hpp"
#include "PointwiseFunctions/AnalyticSolutions/ScalarTensor/Solutions.hpp"
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

namespace ScalarTensor::Solutions {
/*!
 * \brief Set the scalar variables to zero on Minkowski space.
 */
class MinkowskiZeroScalar : /* public evolution::initial_data::InitialData, */
                            /* Why does it work with `virtual`? */
                            public virtual evolution::initial_data::InitialData,
                            public AnalyticSolution,
                            public MarkAsAnalyticSolution {
 public:

  /// The amplitude of the scalar field
  struct Amplitude {
    using type = double;
    static constexpr Options::String help = {
        "The constant pressure throughout the fluid."};
  };

  using options = tmpl::list<Amplitude>;
  static constexpr Options::String help = {
      "Zero scalar field in Minkowski space."};

  MinkowskiZeroScalar() = default;
  MinkowskiZeroScalar(const MinkowskiZeroScalar& /*rhs*/) = default;
  MinkowskiZeroScalar& operator=(const MinkowskiZeroScalar& /*rhs*/) = default;
  MinkowskiZeroScalar(MinkowskiZeroScalar&& /*rhs*/) = default;
  MinkowskiZeroScalar& operator=(MinkowskiZeroScalar&& /*rhs*/) = default;
  ~MinkowskiZeroScalar() override = default;

  MinkowskiZeroScalar(double amplitude);

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit MinkowskiZeroScalar(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(MinkowskiZeroScalar);
  /// \endcond

  // The extra tags below can also be added as common to all scalar tensor
  // solutions
  template <typename DataType, typename Frame = Frame::Inertial>
  using tags = tmpl::flatten<tmpl::list<
      typename AnalyticSolution::template tags<DataType>,
      gr::Tags::DerivDetSpatialMetric<DataType, 3_st, Frame>,
      gr::Tags::TraceExtrinsicCurvature<DataType>,
      gr::Tags::SpatialChristoffelFirstKind<DataType, 3_st, Frame>,
      gr::Tags::SpatialChristoffelSecondKind<DataType, 3_st, Frame>,
      gr::Tags::TraceSpatialChristoffelSecondKind<DataType, 3_st, Frame>>>;

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
  friend bool operator==(const MinkowskiZeroScalar& lhs,
                         const MinkowskiZeroScalar& rhs);

  double amplitude_ = std::numeric_limits<double>::signaling_NaN();
  gr::Solutions::Minkowski<3> background_spacetime_{};
};

bool operator!=(const MinkowskiZeroScalar& lhs, const MinkowskiZeroScalar& rhs);

}  // namespace ScalarTensor::Solutions
