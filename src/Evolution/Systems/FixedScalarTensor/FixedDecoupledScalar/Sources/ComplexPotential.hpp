// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#pragma once

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/CurvedScalarWave/Tags.hpp"
#include "Evolution/Systems/FixedScalarTensor/ScalarDriver/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/ScalarSource.hpp"
#include "Evolution/Systems/ScalarTensor/Sources/Tags.hpp"
#include "Evolution/Systems/ScalarTensor/Tags.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylElectric.hpp"
#include "PointwiseFunctions/GeneralRelativity/WeylMagnetic.hpp"
#include "Utilities/Gsl.hpp"

namespace fe::DecoupledScalar {

void compute_potential_re_part(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    // Real part
    const Scalar<DataVector>& psi_re,
    // Imaginary part
    const Scalar<DataVector>& psi_im, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi) {
  scalar_source->get() = square(mass_psi) * get(psi_re);
}

void compute_potential_im_part(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    // Real part
    const Scalar<DataVector>& psi_re,
    // Imaginary part
    const Scalar<DataVector>& psi_im, const double first_coupling_psi,
    const double second_coupling_psi, const double mass_psi) {
  scalar_source->get() = square(mass_psi) * get(psi_im);
}

// V = (lambda / 4)(Psi_bar * Psi - v^2/ 2)^2 + m^2 Psi_bar Psi
void compute_higgs_potential_re_part(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    // Real part
    const Scalar<DataVector>& psi_re,
    // Imaginary part
    const Scalar<DataVector>& psi_im, const double lambda_parameter,
    const double vev_parameter, const double mass_psi) {
  scalar_source->get() = 0.5 * lambda_parameter *
                         (square(get(psi_re)) + square(get(psi_im))) *
                         get(psi_re);
  scalar_source->get() +=
      -0.5 * lambda_parameter * vev_parameter * vev_parameter * get(psi_re);
  scalar_source->get() += square(mass_psi) * get(psi_re);
}

void compute_higgs_potential_im_part(
    const gsl::not_null<Scalar<DataVector>*> scalar_source,
    // Real part
    const Scalar<DataVector>& psi_re,
    // Imaginary part
    const Scalar<DataVector>& psi_im, const double lambda_parameter,
    const double vev_parameter, const double mass_psi) {
  scalar_source->get() = 0.5 * lambda_parameter *
                         (square(get(psi_re)) + square(get(psi_im))) *
                         get(psi_im);
  scalar_source->get() +=
      -0.5 * lambda_parameter * vev_parameter * vev_parameter * get(psi_im);
  scalar_source->get() += square(mass_psi) * get(psi_im);
}

namespace Tags {

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <typename Frame, typename DataType>
struct ReSourceCompute : ReSource, db::ComputeTag {
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Psi, fe::ScalarDriver::Tags::Psi,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter,
                 ScalarTensor::Tags::ScalarMass>;
  using return_type = Scalar<DataType>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> scalar_source,
      // Real part
      const Scalar<DataVector>& psi_re,
      // Imaginary part
      const Scalar<DataVector>& psi_im, const double first_coupling_psi,
      const double second_coupling_psi,
      const double mass_psi) = &compute_higgs_potential_re_part;
  using base = ReSource;
};

/*!
 * \brief Compute tag for the scalar driver source.
 *
 * \details Call compute_scalar_driver_source.
 */
template <typename Frame, typename DataType>
struct ImSourceCompute : ImSource, db::ComputeTag {
  using argument_tags =
      tmpl::list<CurvedScalarWave::Tags::Psi, fe::ScalarDriver::Tags::Psi,
                 ScalarTensor::Tags::ScalarFirstCouplingParameter,
                 ScalarTensor::Tags::ScalarSecondCouplingParameter,
                 ScalarTensor::Tags::ScalarMass>;
  using return_type = Scalar<DataType>;
  static constexpr void (*function)(
      const gsl::not_null<Scalar<DataVector>*> scalar_source,
      // Real part
      const Scalar<DataVector>& psi_re,
      // Imaginary part
      const Scalar<DataVector>& psi_im, const double first_coupling_psi,
      const double second_coupling_psi,
      const double mass_psi) = &compute_higgs_potential_im_part;
  using base = ImSource;
};

}  // namespace Tags
}  // namespace fe::DecoupledScalar

namespace ScalarTensor::Tags {

struct ReSourceMirrorCompute : ScalarSource, db::ComputeTag {
  using argument_tags = tmpl::list<fe::DecoupledScalar::Tags::ReSource>;
  using base = ScalarSource;
  using return_type = typename base::type;
  static constexpr void function(const gsl::not_null<return_type*> result,
                                 const return_type& csw_var) {
    for (size_t i = 0; i < csw_var.size(); ++i) {
      make_const_view(make_not_null(&std::as_const((*result)[i])), csw_var[i],
                      0, csw_var[i].size());
    }
  }
};
}  // namespace ScalarTensor::Tags

namespace fe::ScalarDriver::Tags {

struct ImSourceMirrorCompute : ScalarDriverSource, db::ComputeTag {
  using argument_tags = tmpl::list<fe::DecoupledScalar::Tags::ImSource>;
  using base = ScalarDriverSource;
  using return_type = typename base::type;
  static constexpr void function(const gsl::not_null<return_type*> result,
                                 const return_type& csw_var) {
    for (size_t i = 0; i < csw_var.size(); ++i) {
      make_const_view(make_not_null(&std::as_const((*result)[i])), csw_var[i],
                      0, csw_var[i].size());
    }
  }
};
}  // namespace fe::ScalarDriver::Tags
