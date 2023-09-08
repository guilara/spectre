// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines DataBox tags for the scalar driver system

#pragma once

#include <cstddef>
#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Options/String.hpp"

/// \cond
class DataVector;
template <class>
class Variables;
/// \endcond

namespace fe::ScalarDriver::Tags {

/*!
 * \brief The scalar field.
 */
struct Psi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "Psi(ScalarDriver)"; }
};

/*!
 * \brief The conjugate momentum of the scalar field.
 *
 * \details Its definition comes from requiring it to be the future-directed
 * time derivative of the scalar field \f$\Psi\f$ in curved spacetime, see
 * \cite Scheel2003vs , Eq. 2.16:
 *
 * \f{align}
 * \Pi :=& -n^a \partial_a \Psi \\
 *     =&  \frac{1}{\alpha}\left(\beta^k \partial_k \Psi -
 * {\partial_t\Psi}\right),\\ \f}
 *
 * where \f$n^a\f$ is the unit normal to spatial slices of the spacetime
 * foliation, \f$\alpha\f$ is the lapse and \f$\beta^i\f$ is the shift vector.
 */
struct Pi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "Pi(ScalarDriver)"; }
};

/*!
 * \brief Auxiliary variable which is analytically the spatial derivative of the
 * scalar field.
 * \details If \f$\Psi\f$ is the scalar field then we define
 * \f$\Phi_{i} = \partial_i \Psi\f$
 */
template<size_t Dim>
struct Phi : db::SimpleTag {
  using type = tnsr::i<DataVector, Dim>;
  static std::string name() { return "Phi(ScalarDriver)"; }
};

struct ConstraintGamma1 : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "ConstraintGamma1(ScalarDriver)"; }
};

struct ConstraintGamma2 : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "ConstraintGamma2(ScalarDriver)"; }
};

/*!
 * \brief Tag for the one-index constraint of the scalar wave
 * system in curved spacetime.
 *
 * For details on how this is defined and computed, see
 * `OneIndexConstraintCompute`.
 */
struct OneIndexConstraint : db::SimpleTag {
  using type = tnsr::i<DataVector, 3_st, Frame::Inertial>;
  static std::string name() { return "OneIndexConstraint(ScalarDriver)"; }
};
/*!
 * \brief Tag for the two-index constraint of the scalar wave
 * system in curved spacetime.
 *
 * For details on how this is defined and computed, see
 * `TwoIndexConstraintCompute`.
 */
struct TwoIndexConstraint : db::SimpleTag {
  using type = tnsr::ij<DataVector, 3_st, Frame::Inertial>;
  static std::string name() { return "TwoIndexConstraint(ScalarDriver)"; }
};

/// @{
/// \brief Tags corresponding to the characteristic fields of the
/// scalar-wave system in curved spacetime.
///
/// \details For details on how these are defined and computed, \see
/// CharacteristicSpeedsCompute
struct VPsi : db::SimpleTag {
  using type = Scalar<DataVector>;
};

template<size_t Dim>
struct VZero : db::SimpleTag {
  using type = tnsr::i<DataVector, Dim, Frame::Inertial>;
};
struct VPlus : db::SimpleTag {
  using type = Scalar<DataVector>;
};
struct VMinus : db::SimpleTag {
  using type = Scalar<DataVector>;
};
/// @}

struct CharacteristicSpeeds : db::SimpleTag {
  using type = std::array<DataVector, 4>;
};

struct LargestCharacteristicSpeed : db::SimpleTag {
  using type = double;
};

struct CharacteristicFields : db::SimpleTag {
  using type = Variables<tmpl::list<VPsi, VZero<3_st>, VPlus, VMinus>>;
};

struct EvolvedFieldsFromCharacteristicFields : db::SimpleTag {
  using type = Variables<tmpl::list<Psi, Pi, Phi<3_st>>>;
};

struct ScalarDriverSource : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct TrackingDiagnostic : db::SimpleTag {
  using type = Scalar<DataVector>;
};

}  // namespace fe::ScalarDriver::Tags

namespace fe::ScalarDriver::OptionTags {
/*!
 * \brief Scalar sigma parameter.
 */
struct ScalarSigmaParameter {
  static std::string name() { return "ScalarDriverSigmaParameter"; }
  using type = double;
  static constexpr Options::String help{
      "Sigma parameter for the scalar diver in code units"};
};

/*!
 * \brief Scalar tau parameter.
 */
struct ScalarTauParameter {
  static std::string name() { return "ScalarDriverTauParameter"; }
  using type = double;
  static constexpr Options::String help{
      "Tau parameter for the scalar driver in code units"};
};

/*!
 * \brief Limiter parameter for the scalar driver source.
 */
struct DriverLimiterParameter {
  static std::string name() { return "DriverLimiterParameter"; }
  using type = double;
  static constexpr Options::String help{
      "Limiter parameter for the scalar driver source"};
  static type lower_bound() { return 0.; }
};

/*!
 * \brief Amplitude of the gaussian function for the constraint damping
 * parameter for the scalar driver equation.
 */
struct AmplitudeConstraintGamma2 {
  static std::string name() {
    return "AmplitudeConstraintGamma2(ScalarDriver)";
  }
  using type = double;
  static constexpr Options::String help{
      "Amplitude of the gaussian function for the constraint damping parameter "
      "for the scalar driver equation."};
  static type lower_bound() { return 0.; }
};

/*!
 * \brief Width of the gaussian function for the constraint damping
 * parameter for the scalar driver equation.
 */
struct SigmaConstraintGamma2 {
  static std::string name() { return "SigmaConstraintGamma2(ScalarDriver)"; }
  using type = double;
  static constexpr Options::String help{
      "Width of the gaussian function for the constraint damping parameter "
      "for the scalar driver equation."};
  static type lower_bound() { return 0.; }
};

/*!
 * \brief Asymptotic value of the gaussian function for the constraint damping
 * parameter for the scalar driver equation.
 */
struct OffsetConstraintGamma2 {
  static std::string name() { return "OffsetConstraintGamma2(ScalarDriver)"; }
  using type = double;
  static constexpr Options::String help{
      "Asymptotic value for the gaussian function for the constraint damping "
      "parameter for the scalar driver equation."};
  static type lower_bound() { return 0.; }
};
}  // namespace fe::ScalarDriver::OptionTags

namespace fe::ScalarDriver::Tags {

/*!
 * \brief The target for the scalar driver field.
 */
struct TargetPsi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "TargetPsi(ScalarDriver)"; }
};

// We should updgrade these to be scalars so that we can prescribe space
// dependence
struct ScalarSigmaParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarSigmaParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double scalar_sigma_parameter) {
    return scalar_sigma_parameter;
  }
};

struct ScalarTauParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarTauParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double scalar_tau_parameter) {
    return scalar_tau_parameter;
  }
};

struct DriverLimiterParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::DriverLimiterParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double limiter_parameter) {
    return limiter_parameter;
  }
};

struct AmplitudeConstraintGamma2 : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::AmplitudeConstraintGamma2>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double amplitude_gaussian) {
    return amplitude_gaussian;
  }
};

struct SigmaConstraintGamma2 : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::SigmaConstraintGamma2>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double sigma_gaussian) {
    return sigma_gaussian;
  }
};

struct OffsetConstraintGamma2 : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::OffsetConstraintGamma2>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double offset_gaussian) {
    return offset_gaussian;
  }
};

}  // namespace fe::ScalarDriver::Tags
