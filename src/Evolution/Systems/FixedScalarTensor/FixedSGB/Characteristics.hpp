// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Characteristics.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarTensorDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Characteristics.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace fe::sgb::Tags {
struct LargestCharacteristicSpeed : db::SimpleTag {
  using type = double;
};

/*!
 * \brief Computes the largest magnitude of the characteristic speeds.
 *
 * \warning Assumes \f$-1\le\gamma_1\le0\f$.
 */
template <typename Frame = Frame::Inertial>
struct ComputeLargestCharacteristicSpeed : db::ComputeTag,
                                           LargestCharacteristicSpeed {
  using argument_tags =
      tmpl::list<::gh::ConstraintDamping::Tags::ConstraintGamma1,
                 gr::Tags::Lapse<DataVector>,
                 gr::Tags::Shift<DataVector, 3, Frame>,
                 gr::Tags::SpatialMetric<DataVector, 3, Frame>,
                 CurvedScalarWave::Tags::ConstraintGamma1>;
  using return_type = double;
  using base = LargestCharacteristicSpeed;
  static void function(const gsl::not_null<double*> speed,
                       // GH arguments
                       const Scalar<DataVector>& gamma_1,
                       const Scalar<DataVector>& lapse,
                       const tnsr::I<DataVector, 3, Frame>& shift,
                       const tnsr::ii<DataVector, 3, Frame>& spatial_metric,
                       // Scalar arguments
                       const Scalar<DataVector>& gamma_1_scalar
                       // Driver arguments
  ) {
    // Largest speed in for ScalarTensor
    double st_largest_speed = 0.0;
    ScalarTensor::Tags::ComputeLargestCharacteristicSpeed<Frame>::function(
        make_not_null(&st_largest_speed), gamma_1, lapse, shift, spatial_metric,
        gamma_1_scalar);
    // Largest speed for ScalarDriver
    double driver_largest_speed = 0.0;
    fe::ScalarTensorDriver::Tags::ComputeLargestCharacteristicSpeed::function(
        make_not_null(&driver_largest_speed), lapse, shift, spatial_metric);
    // Compute the maximum speed
    *speed = std::max(st_largest_speed, driver_largest_speed);
  }
};
}  // namespace fe::sgb::Tags
