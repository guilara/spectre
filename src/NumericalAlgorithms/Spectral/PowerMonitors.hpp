// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

// Implement a function that computes the power monitors
// Take as input a DataVector and a Mesh
// Return std::array<DataVector, Dim> withe the power monitors
// in each dimension

// Overload with input Variables

// Expose the function to Python for use in post-processing scripts

// Add a compute item.
// Add it to the DataBox to compute the power monitors during a simulation
// Output the power monitors to H5 files as 1D volume data
// Add Python script to read the H5 files and plot the monitors
