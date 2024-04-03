#!/usr/bin/env python

# Distributed under the MIT License.
# See LICENSE.txt for details.

import logging
import os
import sys

import click
import h5py
import matplotlib.pyplot as plt
import numpy as np
import rich
import scipy
from scipy import io, optimize

from spectre.Visualization.Plot import (
    apply_stylesheet_command,
    show_or_save_plot_command,
)
from spectre.Visualization.ReadH5 import available_subfiles

logger = logging.getLogger(__name__)


# Read .dat files from Reductions data: Inertial centers, Mass, Spin
def extract_data_from_file(h5_file, subfile_name, functions, x_axis):
    """Extract data from '.dat' datasets in H5 files"""

    with h5py.File(h5_file, "r") as h5file:
        # Print available subfiles and exit
        if subfile_name is None:
            import rich.columns

            rich.print(
                rich.columns.Columns(
                    available_subfiles(h5file, extension=".dat")
                )
            )
            return

        # Open subfile
        if not subfile_name.endswith(".dat"):
            subfile_name += ".dat"
        dat_file = h5file.get(subfile_name)
        if dat_file is None:
            raise click.UsageError(
                f"Unable to open dat file '{subfile_name}'. Available "
                f"files are:\n {available_subfiles(h5file, extension='.dat')}"
            )

        # Read legend from subfile
        legend = list(dat_file.attrs["Legend"])

        # Select x-axis
        if x_axis is None:
            x_axis = legend[0]
        elif x_axis not in legend:
            raise click.UsageError(
                f"Unknown x-axis '{x_axis}'. Available columns are: {legend}"
            )

        num_obs = len(dat_file[:, legend.index(x_axis)])
        out_table = np.reshape(dat_file[:, legend.index(x_axis)], (num_obs, 1))

        # Assemble in table
        for function, label in functions:
            if function not in legend:
                raise click.UsageError(
                    f"Unknown function '{function}'. "
                    f"Available functions are: {legend}"
                )

            out_table = np.append(
                out_table,
                np.reshape(dat_file[:, legend.index(function)], (num_obs, 1)),
                axis=1,
            )

        return out_table


# Compute separation norm
def compute_separation(h5_file, subfile_name_aha, subfile_name_ahb):
    """Compute coordinate separation"""

    functions = [
        ["InertialCenter_x", "x"],
        ["InertialCenter_y", "y"],
        ["InertialCenter_z", "z"],
    ]
    x_axis = "Time"

    # Extract data
    ObjectA_centers = extract_data_from_file(
        h5_file=h5_file,
        subfile_name=subfile_name_aha,
        functions=functions,
        x_axis=x_axis,
    )

    ObjectB_centers = extract_data_from_file(
        h5_file=h5_file,
        subfile_name=subfile_name_ahb,
        functions=functions,
        x_axis=x_axis,
    )

    if (
        subfile_name_aha is None
        or subfile_name_ahb is None
        or subfile_name_aha == subfile_name_ahb
    ):
        raise click.UsageError(
            f"Dat files '{subfile_name_aha}' and '{subfile_name_aha}' are the"
            " same or at least one of them is missing. Choose files for"
            " different objects. Available files are:\n"
            f" {available_subfiles(h5_file, extension='.dat')}"
        )

    # Compute separation
    num_obs = len(ObjectA_centers[:, 0])

    if len(ObjectB_centers[:, 0]) < num_obs:
        num_obs = len(ObjectB_centers[:, 0])

    # Separation vector
    separation_vec = (
        ObjectA_centers[:num_obs, 1:] - ObjectB_centers[:num_obs, 1:]
    )

    # Compute separation norm
    separation_norm = np.zeros((num_obs, 2))
    separation_norm[:, 0] = ObjectA_centers[:, 0]
    separation_norm[:, 1] = np.linalg.norm(separation_vec, axis=1)

    return separation_norm


def windowed_time_derivative_of_separation(data, tmin=None, tmax=None):
    traw = data[:, 0]
    sraw = data[:, 1]

    # Compute separation derivative
    dsdtraw = (sraw[2:] - sraw[0:-2]) / (traw[2:] - traw[0:-2])

    trawcut = traw[1:-1]

    # Apply time window: select values in [tmin, tmax]
    if tmin == None and tmax == None:
        if traw[-1] < 200:
            which_indices = trawcut > 20
        else:
            which_indices = trawcut > 60
    elif tmax == None:
        which_indices = trawcut > tmin
    else:
        which_indices = np.logical_and(trawcut > tmin, trawcut < tmax)

    dsdt = dsdtraw[which_indices]
    t = trawcut[which_indices]

    return t, dsdt


def fit(x, y, model):
    """Fit coordinate separation"""
    F = model["function"]
    inparams = model["initial guess"]

    errfunc = lambda p, x, y: F(p, x) - y
    p, success = optimize.leastsq(errfunc, inparams[:], args=(x, y))

    # Compute rms error of fit
    e2 = (errfunc(p, x, y)) ** 2
    rms = np.sqrt(sum(e2) / np.size(e2))

    return dict([("parameters", p), ("rms", rms), ("success", success)])


def eccentricity_control_updates(
    x, y, model, initial_separation, initial_xcts_values=None
):
    """Compute updates for eccentricity control"""
    fit_results = fit(x=x, y=y, model=model)

    amplitude, omega, phase = fit_results["parameters"][:3]

    # Compute updates for Omega and expansion and compute eccentricity
    dOmg = amplitude / 2.0 / initial_separation * np.sin(phase)
    dadot = -amplitude / initial_separation * np.cos(phase)
    ecc = amplitude / initial_separation / omega

    fit_results["eccentricity"] = ecc
    fit_results["xcts updates"] = dict(
        [("omega update", dOmg), ("expansion update", dadot)]
    )

    # Update xcts parameters if given
    if initial_xcts_values is not None:
        xcts_omega, xcts_expansion = initial_xcts_values
        fit_results["previous xcts values"] = dict(
            [("omega", xcts_omega), ("expansion", xcts_expansion)]
        )
        updated_xcts_omega = xcts_omega + dOmg
        updated_xcts_expansion = xcts_expansion + dadot
        fit_results["updated xcts values"] = dict(
            [
                ("omega", updated_xcts_omega),
                ("expansion", updated_xcts_expansion),
            ]
        )

    return fit_results


def eccentricity_control_digest(x, functions, output=None):
    """Plot ooutput for eccentricity control"""

    for func in functions.items():
        name = func["label"]
        rms = func["fit result"]["rms"]

        logger.info(
            f"==== Function fitted to dOmega/dt: {name:30s},  rms = {rms:4.3g} "
            " ===="
        )

        style = "--"
        F = func["function"]
        ecc = func["fit result"]["eccentricity"]

        errfunc = lambda p, x, y: F(p, x) - y

        if output is not None:
            # Plot dD/dt
            plt.subplot(2, 2, 1)
            plt.plot(
                x,
                F(p, x),
                style,
                label=f"{name:s} \n rms = {rms:2.1e}, ecc = {ecc:4.5f}",
            )
            plt.legend(loc=(1.1, -1.3))

            # Plot residual
            plt.subplot(2, 2, 3)
            plt.plot(x, errfunc(p, x, y), style, label=name)
            plt.title("Residual")

        # Print fit parameters
        p = func["fit result"]["parameters"]
        logger.info("Fit parameters:")
        if np.size(p) == 3:
            logger.info(
                f"Oscillatory part: (B, w)=({p[0]:4.3g}, {p[1]:7.4f}), ",
                f"Polynomial part: ({p[2]:4.2g})",
            )
        else:
            logger.info(
                f"Oscillatory part: (B, w, phi)=({p[0]:4.3g}, {p[1]:6.4f},"
                f" {p[2]:6.4f}), "
            )
            if np.size(p) >= 4:
                logger.info(f"Polynomial part: ({p[3]:4.2g}, ")
                for q in p[4:]:
                    logger.info(f"{q:4.2g}")
                logger.info(")")

        # Print eccentricity
        logger.info(f"Eccentricity based on fit: {ecc:9.6f}")
        # Print suggested updates based on fit
        xcts_updates = func["fit result"]["xcts updates"]
        dOmg = xcts_updates["omega update"]
        dadot = xcts_updates["expansion update"]
        logger.info("Suggested updates based on fit:")
        logger.info(f"(dOmega, dadot) = ({dOmg:+13.10f}, {dadot:+8.6g})")

        if "updated xcts values" in func:
            xcts_omega = func["updated xcts values"]["omega"]
            xcts_expansion = func["updated xcts values"]["expansion"]
            logger.info("Updated Xcts values based on fit:")
            logger.info(
                f"(Omega, adot) = ({(xcts_omega + dOmg):13.10f},"
                f" {(xcts_expansion + dadot):13.10g})"
            )

    return


def eccentricity_control(
    h5_file,
    subfile_name_aha,
    subfile_name_ahb,
    tmin,
    tmax,
    angular_velocity_from_xcts,
    expansion_from_xcts,
    output=None,
):
    """Compute updates based on fits to the coordinate separation for manual
    eccentricity control"""

    data = compute_separation(
        h5_file=h5_file,
        subfile_name_aha=subfile_name_aha,
        subfile_name_ahb=subfile_name_ahb,
    )

    # Separation data (unwindowed)
    sraw = data[:, 1]

    # Compute derivative in time window
    t, dsdt = windowed_time_derivative_of_separation(
        data=data, tmin=tmin, tmax=tmax
    )

    # Collect initial xcts values
    if (
        angular_velocity_from_xcts is not None
        and expansion_from_xcts is not None
    ):
        initial_xcts_values = (
            angular_velocity_from_xcts,
            expansion_from_xcts,
        )
    else:
        initial_xcts_values = None

    # Define functions to fit
    functions = dict([])

    # ==== Restricted fit ====
    functions["F1"] = dict(
        [
            ("label", "B*cos(w*t+np.pi/2)+const"),
            (
                "function",
                lambda p, t: p[0] * np.cos(p[1] * t + np.pi / 2) + p[3],
            ),
            ("initial guess", [0, 0.010, np.pi / 2, 0]),
        ]
    )

    # ==== const + cos ====
    functions["F2"] = dict(
        [
            ("label", "B*cos(w*t+phi)+const"),
            ("function", lambda p, t: p[0] * np.cos(p[1] * t + p[2]) + p[3]),
            ("initial guess", [0, 0.010, 0, 0]),
        ]
    )

    # ==== linear + cos ====
    functions["F3"] = dict(
        [
            ("label", "B*cos(w*t+phi)+linear"),
            (
                "function",
                lambda p, t: p[3] + p[4] * t + p[0] * np.cos(p[1] * t + p[2]),
            ),
            (
                "initial guess",
                [
                    0,
                    0.017,
                    0,
                    0,
                    0,
                ],
            ),
        ]
    )

    # ==== quadratic + cos ====
    functions["F4"] = dict(
        [
            ("label", "B*cos(w*t+phi)+quadratic"),
            (
                "function",
                lambda p, t: p[3]
                + p[4] * t
                + p[5] * t**2
                + p[0] * np.cos(p[1] * t + p[2]),
            ),
            (
                "initial guess",
                [
                    0,
                    0.017,
                    0,
                    0,
                    0,
                ],
            ),
        ]
    )

    # ==== Restricted fit ====
    functions["F1"]["fit results"] = eccentricity_control_updates(
        x=t,
        y=dsdt,
        model=functions["F1"],
        initial_separation=sraw[0],
        initial_xcts_values=initial_xcts_values,
    )

    # ==== const + cos ====
    functions["F2"]["fit results"] = eccentricity_control_updates(
        x=t,
        y=dsdt,
        model=functions["F2"],
        initial_separation=sraw[0],
        initial_xcts_values=initial_xcts_values,
    )

    # ==== linear + cos ====
    functions["F3"]["fit results"] = eccentricity_control_updates(
        x=t,
        y=dsdt,
        model=functions["F3"],
        initial_separation=sraw[0],
        initial_xcts_values=initial_xcts_values,
    )

    # ==== quadratic + cos ====
    # Replace the initial guess with that of the previous solve
    iguess_len = len(functions["F4"]["initial guess"])
    functions["F4"]["initial guess"] = functions["F3"]["fit results"][
        "parameters"
    ][0 : len(iguess_len)]

    functions["F4"]["fit_results"] = eccentricity_control_updates(
        x=t,
        y=dsdt,
        model=functions["F4"],
        initial_separation=sraw[0],
        initial_xcts_values=initial_xcts_values,
    )

    # Print results and plot
    eccentricity_control_digest(x=t, functions=functions, output=output)

    return functions


@click.command(name="eccentricity-control")
@click.argument(
    "h5_file",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
)
@click.option(
    "--subfile-name-aha",
    "-A",
    help=(
        "Name of subfile containing the apparent horizon centers for object A."
    ),
)
@click.option(
    "--subfile-name-ahb",
    "-B",
    help=(
        "Name of subfile containing the apparent horizon centers for object B."
    ),
)
@click.option(
    "--tmin",
    type=float,
    nargs=1,
    help=(
        "The lower time bound to start the fit. Used to remove initial junk and"
        "transients in the coordinate separations."
    ),
)
@click.option(
    "--tmax",
    type=float,
    nargs=1,
    help=(
        "The upper time bound to start the fit. A reasonable value would"
        " include 2-3 orbits."
    ),
)
@click.option(
    "--angular-velocity-from-xcts",
    type=float,
    nargs=1,
    help="Value of the angular velocity used in the Xcts file.",
)
@click.option(
    "--expansion-from-xcts",
    type=float,
    nargs=1,
    help="Value of the expansion velocity (adot) used in the Xcts file.",
)
@apply_stylesheet_command()
@show_or_save_plot_command()
def eccentricity_control_command(
    h5_file,
    subfile_name_aha,
    subfile_name_ahb,
    tmin,
    tmax,
    angular_velocity_from_xcts,
    expansion_from_xcts,
    output,
):
    """Compute updates based on fits to the coordinate separation for manual
    eccentricity control

    Usage:

    Select an appropriate time window without large initial transients and about
    2 to 3 orbits of data. This script uses the coordinate separations between
    Objects A and B to compute a finite difference approximation to the time
    derivative of the orbital velocity (Omega). It then fits different models
    to it.

    For each model, the suggested updates dOmega and dadot based on Newtonian
    estimates are printed. Note that when all models fit the data adequately,
    their updates are similar. When they differ, examine the output plot to
    find a model that is good fit and has small residuals (especially at early
    times).

    Finally, add the selected dOmega and dadot updates to the angular velocity
    and expansion parameters (respectively) in the Xcts input file.

    See ArXiv:gr-qc/0702106 and ArXiv:0710.0158 for more details.

    Limitations:

    1) These eccentricity updates work only for non-precessing binaries.
    2) The time window is manually specified by the user.
    3) The coordinate separation is used, instead of the proper distance.

    See OmegaDoEccRemoval.py in SpEC for improved eccentricity control.

    """
    functions = eccentricity_control(
        h5_file=h5_file,
        subfile_name_aha=subfile_name_aha,
        subfile_name_ahb=subfile_name_ahb,
        tmin=tmin,
        tmax=tmax,
        angular_velocity_from_xcts=angular_velocity_from_xcts,
        expansion_from_xcts=expansion_from_xcts,
        output=output,
    )

    return
