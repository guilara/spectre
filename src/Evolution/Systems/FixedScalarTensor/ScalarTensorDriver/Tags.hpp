// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines DataBox tags for the scalar tensor driver system

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

namespace fe::ScalarTensorDriver::Tags {

/*!
 * \brief The scalar field.
 */
struct Psi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "Psi(ScalarDriver)"; }
};

/*!
 * \brief The conjugate momentum of the scalar field.
 */
struct PiScalar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "Pi(ScalarDriver)"; }
};

/*!
 * \brief Auxiliary variable which is analytically the spatial derivative of the
 * scalar field.
 */
template <size_t Dim>
struct PhiScalar : db::SimpleTag {
  using type = tnsr::i<DataVector, Dim>;
  static std::string name() { return "Phi(ScalarDriver)"; }
};

/*!
 * \brief Tensor driver.
 */
template <typename DataType, size_t Dim, typename Frame>
struct TensorDriver : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
  static std::string name() { return "TensorDriver"; }
};

/*!
 * \brief Conjugate momentum to the tensor driver.
 */
template <typename DataType, size_t Dim, typename Frame>
struct Pi : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
  static std::string name() { return "Phi(TensorDriver)"; }
};

/*!
 * \brief Auxiliary variable which is analytically the spatial derivative of the
 * tensor driver
 */
template <typename DataType, size_t Dim, typename Frame>
struct Phi : db::SimpleTag {
  using type = tnsr::iaa<DataType, Dim, Frame>;
  static std::string name() { return "Phi(TensorDriver)"; }
};

struct ScalarDriverSource : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct ScalarTrackingDiagnostic : db::SimpleTag {
  using type = Scalar<DataVector>;
};

struct TensorDriverSource : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
};

struct TensorTrackingDiagnostic : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
};

/*!
 * \brief Scalar sigma parameter.
 */
struct SigmaParameter : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "DriverSigmaParameter"; }
};

/*!
 * \brief Scalar tau parameter.
 */
struct TauParameter : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "DriverTauParameter"; }
};

/*!
 * \brief The target for the scalar driver field.
 */
struct TargetPsi : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "TargetPsi(ScalarDriver)"; }
};

/*!
 * \brief The target for the scalar driver field.
 */
struct TargetTensorDriver : db::SimpleTag {
  using type = tnsr::aa<DataVector, 3>;
  static std::string name() { return "Target(TensorDriver)"; }
};

}  // namespace fe::ScalarTensorDriver::Tags
