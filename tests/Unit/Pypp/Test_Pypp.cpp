// Distributed under the MIT License.
// See LICENSE.txt for detai

#include <array>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "tests/Unit/Pypp/CheckWithRandomValues.hpp"
#include "tests/Unit/Pypp/Pypp.hpp"
#include "tests/Unit/Pypp/PyppFundamentals.hpp"
#include "tests/Unit/Pypp/SetupLocalPythonEnvironment.hpp"
#include "tests/Unit/TestingFramework.hpp"
#include "tests/Utilities/MakeWithRandomValues.hpp"

// IWYU pragma: no_forward_declare Tensor

SPECTRE_TEST_CASE("Unit.Pypp.none", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  pypp::call<pypp::None>("PyppPyTests", "test_none");
  CHECK_THROWS(pypp::call<pypp::None>("PyppPyTests", "test_numeric", 1, 2));
  CHECK_THROWS(pypp::call<std::string>("PyppPyTests", "test_none"));
}

SPECTRE_TEST_CASE("Unit.Pypp.std::string", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const auto ret = pypp::call<std::string>("PyppPyTests", "test_string",
                                           std::string("test string"));
  CHECK(ret == std::string("back test string"));
  CHECK_THROWS(pypp::call<std::string>("PyppPyTests", "test_string",
                                       std::string("test string_")));
  CHECK_THROWS(pypp::call<double>("PyppPyTests", "test_string",
                                  std::string("test string")));
}

SPECTRE_TEST_CASE("Unit.Pypp.int", "[Pypp][Unit]") {
  /// [pypp_int_test]
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const auto ret = pypp::call<long>("PyppPyTests", "test_numeric", 3, 4);
  CHECK(ret == 3 * 4);
  /// [pypp_int_test]
  CHECK_THROWS(pypp::call<double>("PyppPyTests", "test_numeric", 3, 4));
  CHECK_THROWS(pypp::call<void*>("PyppPyTests", "test_numeric", 3, 4));
  CHECK_THROWS(pypp::call<long>("PyppPyTests", "test_none"));
}

SPECTRE_TEST_CASE("Unit.Pypp.long", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const auto ret = pypp::call<long>("PyppPyTests", "test_numeric", 3l, 4l);
  CHECK(ret == 3l * 4l);
  CHECK_THROWS(pypp::call<double>("PyppPyTests", "test_numeric", 3l, 4l));
  CHECK_THROWS(pypp::call<long>("PyppPyTests", "test_numeric", 3.0, 3.74));
}

SPECTRE_TEST_CASE("Unit.Pypp.unsigned_long", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const auto ret =
      pypp::call<unsigned long>("PyppPyTests", "test_numeric", 3ul, 4ul);
  CHECK(ret == 3ul * 4ul);
  CHECK_THROWS(pypp::call<double>("PyppPyTests", "test_numeric", 3ul, 4ul));
  CHECK_THROWS(
      pypp::call<unsigned long>("PyppPyTests", "test_numeric", 3.0, 3.74));
}

SPECTRE_TEST_CASE("Unit.Pypp.double", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const auto ret =
      pypp::call<double>("PyppPyTests", "test_numeric", 3.49582, 3);
  CHECK(ret == 3.0 * 3.49582);
  CHECK_THROWS(pypp::call<long>("PyppPyTests", "test_numeric", 3.8, 3.9));
  CHECK_THROWS(pypp::call<double>("PyppPyTests", "test_numeric", 3ul, 3ul));
}

SPECTRE_TEST_CASE("Unit.Pypp.std::vector", "[Pypp][Unit]") {
  /// [pypp_vector_test]
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const auto ret = pypp::call<std::vector<double>>(
      "PyppPyTests", "test_vector", std::vector<double>{1.3, 4.9},
      std::vector<double>{4.2, 6.8});
  CHECK(approx(ret[0]) == 1.3 * 4.2);
  CHECK(approx(ret[1]) == 4.9 * 6.8);
  /// [pypp_vector_test]
  CHECK_THROWS(pypp::call<std::string>("PyppPyTests", "test_vector",
                                       std::vector<double>{1.3, 4.9},
                                       std::vector<double>{4.2, 6.8}));
}

SPECTRE_TEST_CASE("Unit.Pypp.std::array", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  // std::arrays and std::vectors should both convert to lists in python so this
  // test calls the same python function as the vector test
  const auto ret = pypp::call<std::array<double, 2>>(
      "PyppPyTests", "test_vector", std::array<double, 2>{{1.3, 4.9}},
      std::array<double, 2>{{4.2, 6.8}});
  CHECK(approx(ret[0]) == 1.3 * 4.2);
  CHECK(approx(ret[1]) == 4.9 * 6.8);
  CHECK_THROWS(pypp::call<double>("PyppPyTests", "test_vector",
                                  std::array<double, 2>{{1.3, 4.9}},
                                  std::array<double, 2>{{4.2, 6.8}}));
}

SPECTRE_TEST_CASE("Unit.Pypp.DataVector", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const auto ret = pypp::call<DataVector>(
      "numpy", "multiply", DataVector{1.3, 4.9}, DataVector{4.2, 6.8});
  CHECK(approx(ret[0]) == 1.3 * 4.2);
  CHECK(approx(ret[1]) == 4.9 * 6.8);
  CHECK_THROWS(pypp::call<std::string>(
      "numpy", "multiply", DataVector{1.3, 4.9}, DataVector{4.2, 6.8}));
  CHECK_THROWS(pypp::call<DataVector>("PyppPyTests", "two_dim_ndarray"));
  CHECK_THROWS(pypp::call<DataVector>("PyppPyTests", "ndarray_of_floats"));
}

SPECTRE_TEST_CASE("Unit.Pypp.Tensor.Double", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};

  const Scalar<double> scalar{0.8};
  const tnsr::A<double, 3> vector{{{3., 4., 5., 6.}}};
  const auto tnsr_ia = []() {
    tnsr::ia<double, 3> tnsr{};
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        tnsr.get(i, j) = i + 2. * j + 1.;
      }
    }
    return tnsr;
  }();
  const auto tnsr_AA = []() {
    tnsr::AA<double, 3> tnsr{};
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        tnsr.get(i, j) = i + j + 1.;
      }
    }
    return tnsr;
  }();
  const auto tnsr_iaa = []() {
    tnsr::iaa<double, 3> tnsr{};
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        for (size_t k = 0; k < 4; ++k) {
          tnsr.get(i, j, k) = 2. * (k + 1) * (j + 1) + i + 1.;
        }
      }
    }
    return tnsr;
  }();
  const auto tnsr_aia = []() {
    tnsr::aia<double, 3> tnsr{};
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        for (size_t k = 0; k < 4; ++k) {
          tnsr.get(i, j, k) = 2. * (k + 1) * (i + 1) + j + 1.5;
        }
      }
    }
    return tnsr;
  }();
  const auto tnsr_aBcc = []() {
    tnsr::aBcc<double, 3> tnsr{};
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        for (size_t k = 0; k < 4; ++k) {
          for (size_t l = 0; l < 4; ++l) {
            tnsr.get(i, j, k, l) = 3. * i + j + (k + 1) * (l + 1) + 1.;
          }
        }
      }
    }
    return tnsr;
  }();

  // Check converting from Tensor to ndarray
  CHECK(pypp::call<bool>("PyppPyTests", "convert_scalar_successful", scalar));
  CHECK(pypp::call<bool>("PyppPyTests", "convert_vector_successful", vector));
  CHECK(pypp::call<bool>("PyppPyTests", "convert_tnsr_ia_successful", tnsr_ia));
  CHECK(pypp::call<bool>("PyppPyTests", "convert_tnsr_AA_successful", tnsr_AA));
  CHECK(
      pypp::call<bool>("PyppPyTests", "convert_tnsr_iaa_successful", tnsr_iaa));
  CHECK(
      pypp::call<bool>("PyppPyTests", "convert_tnsr_aia_successful", tnsr_aia));
  CHECK(pypp::call<bool>("PyppPyTests", "convert_tnsr_aBcc_successful",
                         tnsr_aBcc));

  // Check converting from ndarray to Tensor
  CHECK(scalar == (pypp::call<Scalar<double>>("PyppPyTests", "scalar")));
  CHECK(vector == (pypp::call<tnsr::A<double, 3>>("PyppPyTests", "vector")));
  CHECK(tnsr_ia == (pypp::call<tnsr::ia<double, 3>>("PyppPyTests", "tnsr_ia")));
  CHECK(tnsr_AA == (pypp::call<tnsr::AA<double, 3>>("PyppPyTests", "tnsr_AA")));
  CHECK(tnsr_iaa ==
        (pypp::call<tnsr::iaa<double, 3>>("PyppPyTests", "tnsr_iaa")));
  CHECK(tnsr_aia ==
        (pypp::call<tnsr::aia<double, 3>>("PyppPyTests", "tnsr_aia")));
  CHECK(tnsr_aBcc ==
        (pypp::call<tnsr::aBcc<double, 3>>("PyppPyTests", "tnsr_aBcc")));

  // Check conversion throws with incorrect rank
  CHECK_THROWS((pypp::call<tnsr::i<double, 3>>("PyppPyTests", "scalar")));
  CHECK_THROWS((pypp::call<tnsr::ij<double, 3>>("PyppPyTests", "vector")));
  CHECK_THROWS((pypp::call<Scalar<double>>("PyppPyTests", "tnsr_AA")));
  // Check conversion throws with correct rank but incorrect dimension
  CHECK_THROWS((pypp::call<tnsr::i<double, 3>>("PyppPyTests", "vector")));
  CHECK_THROWS((pypp::call<tnsr::iaa<double, 3>>("PyppPyTests", "tnsr_aia")));
}

SPECTRE_TEST_CASE("Unit.Pypp.Tensor.DataVector", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const size_t npts = 5;
  const Scalar<DataVector> scalar{DataVector(npts, 0.8)};

  const tnsr::A<DataVector, 3> vector{
      {{DataVector(npts, 3.), DataVector(npts, 4.), DataVector(npts, 5.),
        DataVector(npts, 6.)}}};

  const auto tnsr_ia = []() {
    tnsr::ia<DataVector, 3> tnsr{};
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        tnsr.get(i, j) = DataVector(npts, i + 2. * j + 1.);
      }
    }
    return tnsr;
  }();
  const auto tnsr_AA = []() {
    tnsr::AA<DataVector, 3> tnsr{};
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        tnsr.get(i, j) = DataVector(npts, i + j + 1.);
      }
    }
    return tnsr;
  }();
  const auto tnsr_iaa = []() {
    tnsr::iaa<DataVector, 3> tnsr{};
    for (size_t i = 0; i < 3; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        for (size_t k = 0; k < 4; ++k) {
          tnsr.get(i, j, k) = DataVector(npts, 2. * (k + 1) * (j + 1) + i + 1.);
        }
      }
    }
    return tnsr;
  }();
  const auto tnsr_aia = []() {
    tnsr::aia<DataVector, 3> tnsr{};
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        for (size_t k = 0; k < 4; ++k) {
          tnsr.get(i, j, k) =
              DataVector(npts, 2. * (k + 1) * (i + 1) + j + 1.5);
        }
      }
    }
    return tnsr;
  }();
  const auto tnsr_aBcc = []() {
    tnsr::aBcc<DataVector, 3> tnsr{};
    for (size_t i = 0; i < 4; ++i) {
      for (size_t j = 0; j < 4; ++j) {
        for (size_t k = 0; k < 4; ++k) {
          for (size_t l = 0; l < 4; ++l) {
            tnsr.get(i, j, k, l) =
                DataVector(npts, 3. * i + j + (k + 1) * (l + 1) + 1.);
          }
        }
      }
    }
    return tnsr;
  }();

  CHECK(scalar ==
        (pypp::call<Scalar<DataVector>>("PyppPyTests", "identity", scalar)));
  CHECK(vector == (pypp::call<tnsr::A<DataVector, 3>>("PyppPyTests", "identity",
                                                      vector)));
  CHECK(tnsr_ia == (pypp::call<tnsr::ia<DataVector, 3>>("PyppPyTests",
                                                        "identity", tnsr_ia)));
  CHECK(tnsr_AA == (pypp::call<tnsr::AA<DataVector, 3>>("PyppPyTests",
                                                        "identity", tnsr_AA)));
  CHECK(tnsr_iaa == (pypp::call<tnsr::iaa<DataVector, 3>>(
                        "PyppPyTests", "identity", tnsr_iaa)));
  CHECK(tnsr_aia == (pypp::call<tnsr::aia<DataVector, 3>>(
                        "PyppPyTests", "identity", tnsr_aia)));
  CHECK(tnsr_aBcc == (pypp::call<tnsr::aBcc<DataVector, 3>>(
                         "PyppPyTests", "identity", tnsr_aBcc)));
}

namespace {
template <typename T>
void test_einsum(const T& used_for_size) {
  std::random_device r;
  const auto seed = r();
  std::mt19937 generator(seed);
  INFO("seed" << seed);
  std::uniform_real_distribution<> dist(-10., 10.);

  const auto nn_generator = make_not_null(&generator);
  const auto nn_dist = make_not_null(&dist);

  const auto scalar =
      make_with_random_values<Scalar<T>>(nn_generator, nn_dist, used_for_size);
  const auto vector = make_with_random_values<tnsr::A<T, 3>>(
      nn_generator, nn_dist, used_for_size);
  const auto tnsr_ia = make_with_random_values<tnsr::ia<T, 3>>(
      nn_generator, nn_dist, used_for_size);
  const auto tnsr_AA = make_with_random_values<tnsr::AA<T, 3>>(
      nn_generator, nn_dist, used_for_size);
  const auto tnsr_iaa = make_with_random_values<tnsr::iaa<T, 3>>(
      nn_generator, nn_dist, used_for_size);

  const auto expected = [&scalar, &vector, &tnsr_ia, &tnsr_AA, &tnsr_iaa,
                         &used_for_size]() {
    auto tnsr = make_with_value<tnsr::i<T, 3>>(used_for_size, 0.);
    for (size_t i = 0; i < 3; ++i) {
      for (size_t a = 0; a < 4; ++a) {
        tnsr.get(i) += scalar.get() * vector.get(a) * tnsr_ia.get(i, a);
        for (size_t b = 0; b < 4; ++b) {
          tnsr.get(i) += tnsr_AA.get(a, b) * tnsr_iaa.get(i, a, b);
        }
      }
    }
    return tnsr;
  }();
  /// [einsum_example]
  const auto tensor_from_python = pypp::call<tnsr::i<T, 3>>(
      "PyppPyTests", "test_einsum", scalar, vector, tnsr_ia, tnsr_AA, tnsr_iaa);
  /// [einsum_example]
  CHECK_ITERABLE_APPROX(expected, tensor_from_python);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Pypp.EinSum", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  test_einsum<double>(0.);
  test_einsum<DataVector>(DataVector(5));
}

namespace {
template <typename T>
void check_single_not_null0(const gsl::not_null<T*> result,
                            const T& t0) noexcept {
  *result = t0 + 5.0;
}

template <typename T>
void check_single_not_null0_scalar(const gsl::not_null<T*> result,
                                   const T& t0) noexcept {
  get(*result) = get(t0) + 5.0;
}

template <typename T>
void check_single_not_null1(const gsl::not_null<T*> result, const T& t0,
                            const T& t1) noexcept {
  *result = t0 + t1;
}

template <typename T>
void check_single_not_null2(const gsl::not_null<T*> result, const T& t0,
                            const T& t1) noexcept {
  *result = sqrt(t0) + 1.0 / sqrt(-t1);
}

template <typename T>
void check_single_not_null1_scalar(const gsl::not_null<T*> result, const T& t0,
                                   const T& t1) noexcept {
  get(*result) = get(t0) + get(t1);
}

template <typename T>
void check_single_not_null2_scalar(const gsl::not_null<T*> result, const T& t0,
                                   const T& t1) noexcept {
  get(*result) = sqrt(get(t0)) + 1.0 / sqrt(-get(t1));
}

template <typename T>
void check_double_not_null0(const gsl::not_null<T*> result0,
                            const gsl::not_null<T*> result1,
                            const T& t0) noexcept {
  *result0 = t0 + 5.0;
  *result1 = 2.0 * t0 + 5.0;
}

template <typename T>
void check_double_not_null0_scalar(const gsl::not_null<T*> result0,
                                   const gsl::not_null<T*> result1,
                                   const T& t0) noexcept {
  get(*result0) = get(t0) + 5.0;
  get(*result1) = 2.0 * get(t0) + 5.0;
}

template <typename T>
void check_double_not_null1(const gsl::not_null<T*> result0,
                            const gsl::not_null<T*> result1, const T& t0,
                            const T& t1) noexcept {
  *result0 = t0 + t1;
  *result1 = 2.0 * t0 + t1;
}

template <typename T>
void check_double_not_null1_scalar(const gsl::not_null<T*> result0,
                                   const gsl::not_null<T*> result1, const T& t0,
                                   const T& t1) noexcept {
  get(*result0) = get(t0) + get(t1);
  get(*result1) = 2.0 * get(t0) + get(t1);
}

template <typename T>
void check_double_not_null2(const gsl::not_null<T*> result0,
                            const gsl::not_null<T*> result1, const T& t0,
                            const T& t1) noexcept {
  *result0 = sqrt(t0) + 1.0 / sqrt(-t1);
  *result1 = 2.0 * t0 + t1;
}

template <typename T>
void check_double_not_null2_scalar(const gsl::not_null<T*> result0,
                                   const gsl::not_null<T*> result1, const T& t0,
                                   const T& t1) noexcept {
  get(*result0) = sqrt(get(t0)) + 1.0 / sqrt(-get(t1));
  get(*result1) = 2.0 * get(t0) + get(t1);
}


template <typename T>
T check_by_value0(const T& t0) noexcept {
  return t0 + 5.0;
}

template <typename T>
T check_by_value0_scalar(const T& t0) noexcept {
  return T{get(t0) + 5.0};
}

template <typename T>
T check_by_value1(const T& t0, const T& t1) noexcept {
  return t0 + t1;
}

template <typename T>
T check_by_value1_scalar(const T& t0, const T& t1) noexcept {
  return T{get(t0) + get(t1)};
}

template <typename T>
T check_by_value2(const T& t0, const T& t1) noexcept {
  return sqrt(t0) + 1.0/sqrt(-t1);
}

template <typename T>
T check_by_value2_scalar(const T& t0, const T& t1) noexcept {
  return T{sqrt(get(t0)) + 1.0 / sqrt(-get(t1))};
}

}  // namespace

SPECTRE_TEST_CASE("Unit.Pypp.CheckWithPython", "[Pypp][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{"Pypp/"};
  const DataVector dv(5);
  const double doub(0.);
  const Scalar<double> scalar_double{5.0};
  const Scalar<DataVector> scalar_dv{size_t{5}};
  pypp::check_with_random_values<1>(&check_single_not_null0<double>,
                                    "PyppPyTests", {"check_single_not_null0"},
                                    {{{-10.0, 10.0}}}, doub);
  pypp::check_with_random_values<1>(&check_single_not_null0<DataVector>,
                                    "PyppPyTests", {"check_single_not_null0"},
                                    {{{-10.0, 10.0}}}, dv);
  pypp::check_with_random_values<1>(&check_single_not_null1<double>,
                                    "PyppPyTests", {"check_single_not_null1"},
                                    {{{-10.0, 10.0}}}, doub);
  pypp::check_with_random_values<1>(&check_single_not_null1<DataVector>,
                                    "PyppPyTests", {"check_single_not_null1"},
                                    {{{-10.0, 10.0}}}, dv);
  pypp::check_with_random_values<2>(&check_single_not_null2<double>,
                                    "PyppPyTests", {"check_single_not_null2"},
                                    {{{0.0, 10.0}, {-10.0, 0.0}}}, doub);
  pypp::check_with_random_values<2>(&check_single_not_null2<DataVector>,
                                    "PyppPyTests", {"check_single_not_null2"},
                                    {{{0.0, 10.0}, {-10.0, 0.0}}}, dv);
  pypp::check_with_random_values<1>(
      &check_single_not_null0_scalar<Scalar<double>>, "PyppPyTests",
      {"check_single_not_null0"}, {{{-10.0, 10.0}}}, scalar_double);
  pypp::check_with_random_values<1>(
      &check_single_not_null0_scalar<Scalar<DataVector>>, "PyppPyTests",
      {"check_single_not_null0"}, {{{-10.0, 10.0}}}, scalar_dv);
  pypp::check_with_random_values<1>(
      &check_single_not_null1_scalar<Scalar<double>>, "PyppPyTests",
      {"check_single_not_null1"}, {{{-10.0, 10.0}}}, scalar_double);
  pypp::check_with_random_values<1>(
      &check_single_not_null1_scalar<Scalar<DataVector>>, "PyppPyTests",
      {"check_single_not_null1"}, {{{-10.0, 10.0}}}, scalar_dv);
  pypp::check_with_random_values<2>(
      &check_single_not_null2_scalar<Scalar<double>>, "PyppPyTests",
      {"check_single_not_null2"}, {{{0.0, 10.0}, {-10.0, 0.0}}}, scalar_double);
  pypp::check_with_random_values<2>(
      &check_single_not_null2_scalar<Scalar<DataVector>>, "PyppPyTests",
      {"check_single_not_null2"}, {{{0.0, 10.0}, {-10.0, 0.0}}}, scalar_dv);

  pypp::check_with_random_values<1>(
      &check_double_not_null0<double>, "PyppPyTests",
      {"check_double_not_null0_result0", "check_double_not_null0_result1"},
      {{{-10.0, 10.0}}}, doub);
  pypp::check_with_random_values<1>(
      &check_double_not_null0<DataVector>, "PyppPyTests",
      {"check_double_not_null0_result0", "check_double_not_null0_result1"},
      {{{-10.0, 10.0}}}, dv);
  pypp::check_with_random_values<1>(
      &check_double_not_null1<double>, "PyppPyTests",
      {"check_double_not_null1_result0", "check_double_not_null1_result1"},
      {{{-10.0, 10.0}}}, doub);
  pypp::check_with_random_values<1>(
      &check_double_not_null1<DataVector>, "PyppPyTests",
      {"check_double_not_null1_result0", "check_double_not_null1_result1"},
      {{{-10.0, 10.0}}}, dv);
  pypp::check_with_random_values<2>(
      &check_double_not_null2<double>, "PyppPyTests",
      {"check_double_not_null2_result0", "check_double_not_null2_result1"},
      {{{0.0, 10.0}, {-10.0, 0.0}}}, doub);
  pypp::check_with_random_values<2>(
      &check_double_not_null2<DataVector>, "PyppPyTests",
      {"check_double_not_null2_result0", "check_double_not_null2_result1"},
      {{{0.0, 10.0}, {-10.0, 0.0}}}, dv);
  pypp::check_with_random_values<1>(
      &check_double_not_null0_scalar<Scalar<double>>, "PyppPyTests",
      {"check_double_not_null0_result0", "check_double_not_null0_result1"},
      {{{-10.0, 10.0}}}, scalar_double);
  pypp::check_with_random_values<1>(
      &check_double_not_null0_scalar<Scalar<DataVector>>, "PyppPyTests",
      {"check_double_not_null0_result0", "check_double_not_null0_result1"},
      {{{-10.0, 10.0}}}, scalar_dv);
  pypp::check_with_random_values<1>(
      &check_double_not_null1_scalar<Scalar<double>>, "PyppPyTests",
      {"check_double_not_null1_result0", "check_double_not_null1_result1"},
      {{{-10.0, 10.0}}}, scalar_double);
  pypp::check_with_random_values<1>(
      &check_double_not_null1_scalar<Scalar<DataVector>>, "PyppPyTests",
      {"check_double_not_null1_result0", "check_double_not_null1_result1"},
      {{{-10.0, 10.0}}}, scalar_dv);
  pypp::check_with_random_values<2>(
      &check_double_not_null2_scalar<Scalar<double>>, "PyppPyTests",
      {"check_double_not_null2_result0", "check_double_not_null2_result1"},
      {{{0.0, 10.0}, {-10.0, 0.0}}}, scalar_double);
  pypp::check_with_random_values<2>(
      &check_double_not_null2_scalar<Scalar<DataVector>>, "PyppPyTests",
      {"check_double_not_null2_result0", "check_double_not_null2_result1"},
      {{{0.0, 10.0}, {-10.0, 0.0}}}, scalar_dv);

  pypp::check_with_random_values<1>(&check_by_value0<double>, "PyppPyTests",
                                    "check_by_value0", {{{-10.0, 10.0}}}, doub);
  pypp::check_with_random_values<1>(&check_by_value0<DataVector>, "PyppPyTests",
                                    "check_by_value0", {{{-10.0, 10.0}}}, dv);
  pypp::check_with_random_values<1>(&check_by_value1<double>, "PyppPyTests",
                                    "check_by_value1", {{{-10.0, 10.0}}}, doub);
  pypp::check_with_random_values<1>(&check_by_value1<DataVector>, "PyppPyTests",
                                    "check_by_value1", {{{-10.0, 10.0}}}, dv);
  pypp::check_with_random_values<2>(&check_by_value2<double>, "PyppPyTests",
                                    "check_by_value2",
                                    {{{0.0, 10.0}, {-10.0, 0.0}}}, doub);
  pypp::check_with_random_values<2>(&check_by_value2<DataVector>, "PyppPyTests",
                                    "check_by_value2",
                                    {{{0.0, 10.0}, {-10.0, 0.0}}}, dv);
  pypp::check_with_random_values<1>(&check_by_value0_scalar<Scalar<double>>,
                                    "PyppPyTests", "check_by_value0",
                                    {{{-10.0, 10.0}}}, scalar_double);
  pypp::check_with_random_values<1>(&check_by_value0_scalar<Scalar<DataVector>>,
                                    "PyppPyTests", "check_by_value0",
                                    {{{-10.0, 10.0}}}, scalar_dv);
  pypp::check_with_random_values<1>(&check_by_value1_scalar<Scalar<double>>,
                                    "PyppPyTests", "check_by_value1",
                                    {{{-10.0, 10.0}}}, scalar_double);
  pypp::check_with_random_values<1>(&check_by_value1_scalar<Scalar<DataVector>>,
                                    "PyppPyTests", "check_by_value1",
                                    {{{-10.0, 10.0}}}, scalar_dv);
  pypp::check_with_random_values<2>(
      &check_by_value2_scalar<Scalar<double>>, "PyppPyTests", "check_by_value2",
      {{{0.0, 10.0}, {-10.0, 0.0}}}, scalar_double);
  pypp::check_with_random_values<2>(&check_by_value2_scalar<Scalar<DataVector>>,
                                    "PyppPyTests", "check_by_value2",
                                    {{{0.0, 10.0}, {-10.0, 0.0}}}, scalar_dv);
}
