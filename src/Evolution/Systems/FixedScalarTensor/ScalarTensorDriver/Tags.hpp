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
template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct TensorDriver : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
  static std::string name() { return "TensorDriver"; }
};

/*!
 * \brief Conjugate momentum to the tensor driver.
 */
template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct Pi : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
  static std::string name() { return "Pi(TensorDriver)"; }
};

/*!
 * \brief Auxiliary variable which is analytically the spatial derivative of the
 * tensor driver
 */
template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
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

template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct TensorDriverSource : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
};

template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct TensorTrackingDiagnostic : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
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
struct TargetScalar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "TargetScalar"; }
};

/*!
 * \brief The target for the scalar driver field.
 */
template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct TargetTensor : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
  static std::string name() { return "TargetTensor"; }
};

/// @{
/// \brief Tags corresponding to the characteristic fields.
///
/// \details For advection drivers, the characteristic fields are the same as
/// the evolved fields
template <typename DataType>
struct VScalarDriver : db::SimpleTag {
  using type = Scalar<DataType>;
};
template <typename DataType>
struct VPiScalar : db::SimpleTag {
  using type = Scalar<DataType>;
};
template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct VTensorDriver : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
};
template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct VPi : db::SimpleTag {
  using type = tnsr::aa<DataType, Dim, Frame>;
};
/// @}

template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct CharacteristicSpeeds : db::SimpleTag {
  using type = std::array<DataType, 4>;
};

struct LargestCharacteristicSpeed : db::SimpleTag {
  using type = double;
};

template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct CharacteristicFields : db::SimpleTag {
  using type =
      Variables<tmpl::list<VScalarDriver<DataType>, VPiScalar<DataType>,
                           VTensorDriver<DataType, Dim, Frame>,
                           VPi<DataType, Dim, Frame>>>;
};

template <typename DataType, size_t Dim, typename Frame = Frame::Inertial>
struct EvolvedFieldsFromCharacteristicFields : db::SimpleTag {
  using type =
      Variables<tmpl::list<Psi, PiScalar, TensorDriver<DataType, Dim, Frame>,
                           Pi<DataType, Dim, Frame>>>;
};

/*!
 * \brief Tensor driver.
 */
struct TensorDriverTrace : db::SimpleTag {
  using type = Scalar<DataVector>;
  static std::string name() { return "TensorDriverTrace"; }
};

/*!
 * \brief Tensor driver.
 */
struct TensorDriverSpatialProjection : db::SimpleTag {
  using type = tnsr::ii<DataVector, 3, Frame::Inertial>;
  static std::string name() { return "TensorDriverSpatialProjection"; }
};

}  // namespace fe::ScalarTensorDriver::Tags
