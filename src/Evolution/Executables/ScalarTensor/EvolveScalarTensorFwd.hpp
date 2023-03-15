// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

/// \cond
namespace GeneralizedHarmonic {
template <size_t Dim>
struct System;
namespace Solutions {
template <typename GrSolution>
struct WrappedGr;
}  // namespace Solutions
}  // namespace GeneralizedHarmonic

namespace gr::Solutions {
template <size_t Dim>
class Minkowski;
class KerrSchild;
}  // namespace gr::Solutions

namespace ScalarWave::Solutions {
template <size_t Dim>
class PlaneWave;
}  // namespace ScalarWave::Solutions

template <typename InitialData>
struct EvolutionMetavars;
/// \endcond
