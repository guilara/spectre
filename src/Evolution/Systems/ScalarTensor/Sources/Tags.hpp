// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "Options/String.hpp"

namespace ScalarTensor {
namespace Tags {
// Extra compute tags for debugging
/*!
 * \brief The \f$ f'(\Psi) \f$ term.
 */
struct CouplingFunctionDerivative : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "CouplingFunctionDerivative"; }
};

/*!
 * \brief The \f$ 8 (E^2 - B^2) \f$ term.
 */
struct GBScalar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "GBScalar"; }
};

/*!
 * \brief The GB scalar term with nonvaccuum contributions from the scalar.
 */
struct OrderReducedGBScalar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "OrderReducedGBScalar"; }
};

/*!
 * \brief Normal normal projection of the second covariant derivative of the
 * scalar.
 */
struct nnDDKG : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "nnDDKG"; }
};

/*!
 * \brief Normal spatial projection of the second covariant derivative of the
 * scalar.
 */
struct nsDDKG : db::SimpleTag {
  using type = tnsr::i<DataVector, 3>;
  static std::string name() { return "nsDDKG"; }
};

/*!
 * \brief Spatial spatial projection of the second covariant derivative of the
 * scalar.
 */
struct ssDDKG : db::SimpleTag {
  using type = tnsr::ij<DataVector, 3>;
  static std::string name() { return "ssDDKG"; }
};

/*!
 * \brief Rhs of the Psi equation.
 */
struct RhsPsi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "Rhs(Csw(Psi))"; }
};

/*!
 * \brief Rhs of the Pi equation.
 */
struct RhsPi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "Rhs(Csw(Pi))"; }
};

/*!
 * \brief Rhs of the Pi equation.
 */
struct RhsPhi : db::SimpleTag {
  using type = tnsr::i<DataVector, 3>;
  static std::string name() { return "Rhs(Csw(Phi))"; }
};

/*!
 * \brief Normal normal projection of the order reduced H tensor.
 */
struct OrderReducednnH : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "OrderReducednnH"; }
};

/*!
 * \brief Normal spatial projection of the order reduced H tensor.
 */
struct OrderReducednsH : db::SimpleTag {
  using type = tnsr::i<DataVector, 3>;
  static std::string name() { return "OrderReducednsH"; }
};

/*!
 * \brief Spatial spatial projection of the order reduced H tensor.
 */
struct OrderReducedssH : db::SimpleTag {
  using type = tnsr::ij<DataVector, 3>;
  static std::string name() { return "OrderReducedssH"; }
};

/*!
 * \brief Order reduced H tensor.
 */
struct OrderReducedHTensor : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
  static std::string name() { return "OrderReducedHTensor"; }
};

/*!
 * \brief Order reduced Q tensor.
 */
struct OrderReducedQTensor : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
  static std::string name() { return "OrderReducedQTensor"; }
};

/*!
 * \brief S cross B for normal spatial projection of order reduced H tensor.
 */
struct SCrossB : db::SimpleTag {
  using type = tnsr::i<DataVector, 3>;
  static std::string name() { return "SCrossB"; }
};

/*!
 * \brief j cross B for spatial spatial projection of order reduced H tensor.
 */
struct JCrossB : db::SimpleTag {
  using type = tnsr::ij<DataVector, 3>;
  static std::string name() { return "JCrossB"; }
};

}  // namespace Tags

namespace OptionTags {

/*!
 * \brief Coupling parameter 1.
 */
struct ScalarFirstCouplingParameter {
  static std::string name() { return "EtaCouplingParameter"; }
  using type = double;
  static constexpr Options::String help{
      "First coupling parameter entering in the scalar field source, in code "
      "units"};
};

/*!
 * \brief Coupling parameter 2.
 */
struct ScalarSecondCouplingParameter {
  static std::string name() { return "ZetaCouplingParameter"; }
  using type = double;
  static constexpr Options::String help{
      "Second coupling parameter entering in the scalar field source, in code "
      "units"};
};
}  // namespace OptionTags

namespace Tags {

struct ScalarFirstCouplingParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarFirstCouplingParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double first_coupling_psi) {
    return first_coupling_psi;
  }
};

struct ScalarSecondCouplingParameter : db::SimpleTag {
  using type = double;
  using option_tags = tmpl::list<OptionTags::ScalarSecondCouplingParameter>;
  static constexpr bool pass_metavariables = false;
  static double create_from_options(const double second_coupling_psi) {
    return second_coupling_psi;
  }
};
}  // namespace Tags
}  // namespace ScalarTensor
