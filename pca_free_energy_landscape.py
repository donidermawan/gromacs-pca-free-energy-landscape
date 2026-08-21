#!/usr/bin/env python3

"""
PCA FREE-ENERGY LANDSCAPE ANALYSIS
==================================

Input:
    pc1_pc2_2d.xvg

Expected XVG format:
    PC1    PC2

Outputs:
    PCA_Free_Energy_Landscape.png
    PCA_Free_Energy_Landscape.pdf
    PCA_Free_Energy_Data.csv

Method:
    1. Read PC1/PC2 coordinates from GROMACS XVG
    2. Generate a 2D population histogram
    3. Apply optional light Gaussian smoothing
    4. Normalize the probability distribution
    5. Calculate relative free energy:

           ΔG = -RT ln(P / Pmax)

       Therefore:
           most populated region = 0 kJ/mol

    6. Remove extremely poorly sampled regions
       from visualization
    7. Clip only the displayed free-energy range
    8. Generate a sharp 3D free-energy landscape
    9. Generate a 2D free-energy contour map
   10. Export the original numerical data

"""


# ============================================================
# IMPORTS
# ============================================================

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter

from matplotlib.colors import Normalize


# ============================================================
# USER SETTINGS
# ============================================================

# ------------------------------------------------------------
# INPUT FILE
# ------------------------------------------------------------

INPUT_FILE = "pc1_pc2_2d.xvg"


# ------------------------------------------------------------
# TEMPERATURE
# ------------------------------------------------------------

TEMPERATURE = 300.0


# ------------------------------------------------------------
# HISTOGRAM RESOLUTION
# ------------------------------------------------------------
#
# Higher value:
#     more spatial detail
#
# Lower value:
#     smoother / less detailed landscape
#
# 100 × 100 is a good starting point.
# ------------------------------------------------------------

NBINS = 100


# ------------------------------------------------------------
# GAUSSIAN SMOOTHING
# ------------------------------------------------------------
#
# 0.0 = no smoothing
# 0.5 = very light smoothing
# 1.0 = light smoothing
# 2.0 = moderate smoothing
# 3.0 = strong smoothing
#
# The reference shape is relatively sharp and rugged.
#
# Therefore 0.5 is recommended.
#
# If you want an even sharper landscape:
#
#     SMOOTHING = 0.0
#
# ------------------------------------------------------------

SMOOTHING = 0.5


# ------------------------------------------------------------
# NUMBER OF CONTOUR LEVELS
# ------------------------------------------------------------

N_LEVELS = 20


# ============================================================
# FREE-ENERGY DISPLAY SETTINGS
# ============================================================

# Options:
#
#     "percentile"
#     "fixed"
#
# percentile:
#     automatically determine the upper displayed
#     free-energy limit.
#
# fixed:
#     manually define MAX_FREE_ENERGY.
# ------------------------------------------------------------

DISPLAY_MODE = "percentile"


# ------------------------------------------------------------
# Used when DISPLAY_MODE = "percentile"
#
# 95 means that the upper visual limit is the
# 95th percentile of the sampled free-energy values.
# ------------------------------------------------------------

ENERGY_PERCENTILE = 95


# ------------------------------------------------------------
# Used when DISPLAY_MODE = "fixed"
# ------------------------------------------------------------

MAX_FREE_ENERGY = 10.0


# ============================================================
# PROBABILITY THRESHOLD
# ============================================================

# Bins with probability below:
#
#     Pmax × MIN_PROBABILITY_FRACTION
#
# are considered poorly sampled for visualization.
#
# Example:
#
#     0.001 = 0.1% of Pmax
#
# IMPORTANT:
#
# This affects visualization only.
# The original calculated values remain in the CSV.
# ------------------------------------------------------------

MIN_PROBABILITY_FRACTION = 0.001


# ============================================================
# 3D VIEW SETTINGS
# ============================================================

# Elevation angle of the 3D landscape.
ELEVATION = 28


# Azimuth angle of the 3D landscape.
AZIMUTH = -125


# ============================================================
# OUTPUT FILES
# ============================================================

OUTPUT_PNG = (
    "PCA_Free_Energy_Landscape.png"
)

OUTPUT_PDF = (
    "PCA_Free_Energy_Landscape.pdf"
)

OUTPUT_CSV = (
    "PCA_Free_Energy_Data.csv"
)


# ============================================================
# PHYSICAL CONSTANT
# ============================================================

# Universal gas constant
#
# kJ mol-1 K-1
# ------------------------------------------------------------

R = 0.008314462618


# ============================================================
# FUNCTION: READ GROMACS XVG
# ============================================================

def read_xvg(filename):

    """
    Read PC1 and PC2 values from a GROMACS XVG file.

    Lines beginning with:
        @
        #

    are ignored.

    Returns:
        pc1 : numpy array
        pc2 : numpy array
    """

    pc1 = []
    pc2 = []

    with open(filename, "r") as file:

        for line in file:

            line = line.strip()

            # ------------------------------------------------
            # Ignore empty lines
            # ------------------------------------------------

            if not line:
                continue


            # ------------------------------------------------
            # Ignore GROMACS metadata
            # ------------------------------------------------

            if line.startswith("@"):
                continue

            if line.startswith("#"):
                continue


            # ------------------------------------------------
            # Remove inline comments
            # ------------------------------------------------

            if "#" in line:

                line = line.split("#")[0]


            # ------------------------------------------------
            # Split columns
            # ------------------------------------------------

            parts = line.split()


            if len(parts) < 2:
                continue


            # ------------------------------------------------
            # Convert first two columns to float
            # ------------------------------------------------

            try:

                x = float(parts[0])
                y = float(parts[1])

            except ValueError:

                continue


            pc1.append(x)
            pc2.append(y)


    # ========================================================
    # CHECK DATA
    # ========================================================

    if len(pc1) == 0:

        raise RuntimeError(
            f"No numerical PC1/PC2 data found in "
            f"{filename}"
        )


    # ========================================================
    # CONVERT TO NUMPY
    # ========================================================

    pc1 = np.asarray(
        pc1,
        dtype=float
    )

    pc2 = np.asarray(
        pc2,
        dtype=float
    )


    return pc1, pc2


# ============================================================
# START ANALYSIS
# ============================================================

print()
print("=" * 70)
print(" PCA FREE-ENERGY LANDSCAPE")
print(" GROMACS PC1 / PC2")
print("=" * 70)
print()


# ============================================================
# READ INPUT
# ============================================================

print("Reading:", INPUT_FILE)

pc1, pc2 = read_xvg(
    INPUT_FILE
)


# ============================================================
# REMOVE INVALID VALUES
# ============================================================

valid = (
    np.isfinite(pc1)
    &
    np.isfinite(pc2)
)


if not np.all(valid):

    removed = np.sum(~valid)

    print()
    print(
        f"Warning: removing {removed} "
        f"invalid PC1/PC2 points."
    )

    pc1 = pc1[valid]
    pc2 = pc2[valid]


if len(pc1) == 0:

    raise RuntimeError(
        "No valid PC1/PC2 data remain."
    )


# ============================================================
# BASIC INFORMATION
# ============================================================

print()
print("Successfully loaded PCA data.")

print(
    f"Number of frames: {len(pc1)}"
)

print(
    f"PC1 range: "
    f"{pc1.min():.5f} "
    f"to "
    f"{pc1.max():.5f}"
)

print(
    f"PC2 range: "
    f"{pc2.min():.5f} "
    f"to "
    f"{pc2.max():.5f}"
)


# ============================================================
# 2D HISTOGRAM
# ============================================================

hist, xedges, yedges = np.histogram2d(

    pc1,
    pc2,

    bins=NBINS

)


print()
print(
    f"2D histogram generated: "
    f"{NBINS} × {NBINS} bins"
)


# ============================================================
# GAUSSIAN SMOOTHING
# ============================================================

if SMOOTHING > 0:

    probability = gaussian_filter(

        hist.astype(float),

        sigma=SMOOTHING,

        mode="nearest"

    )

else:

    probability = hist.astype(float)


# ============================================================
# NORMALIZE PROBABILITY
# ============================================================

total_probability = probability.sum()


if total_probability <= 0:

    raise RuntimeError(
        "Probability distribution is empty."
    )


probability /= total_probability


# ============================================================
# MAXIMUM PROBABILITY
# ============================================================

Pmax = probability.max()


if Pmax <= 0:

    raise RuntimeError(
        "Maximum probability is zero."
    )


# ============================================================
# FREE-ENERGY CALCULATION
# ============================================================

"""
Relative free energy:

    ΔG = -RT ln(P/Pmax)

where:

    R    = gas constant
    T    = temperature
    P    = probability
    Pmax = maximum probability

Therefore:

    P = Pmax

gives:

    ΔG = 0 kJ/mol

The most populated conformational region
is therefore the global free-energy minimum.
"""

free_energy = np.full(

    probability.shape,

    np.nan,

    dtype=float

)


positive_probability = (
    probability > 0
)


free_energy[
    positive_probability
] = (

    -R
    * TEMPERATURE
    * np.log(
        probability[
            positive_probability
        ]
        / Pmax
    )

)


# ============================================================
# GRID CENTERS
# ============================================================

xcenters = (

    xedges[:-1]
    + xedges[1:]

) / 2.0


ycenters = (

    yedges[:-1]
    + yedges[1:]

) / 2.0


X, Y = np.meshgrid(

    xcenters,
    ycenters,

    indexing="ij"

)


# ============================================================
# PROBABILITY THRESHOLD
# ============================================================

probability_threshold = (

    Pmax
    * MIN_PROBABILITY_FRACTION

)


sampled_mask = (

    probability
    >= probability_threshold

)


# ============================================================
# MASK POORLY SAMPLED REGIONS
# ============================================================

free_energy_masked = (
    free_energy.copy()
)


free_energy_masked[
    ~sampled_mask
] = np.nan


# ============================================================
# EXTRACT FINITE FREE-ENERGY VALUES
# ============================================================

finite_energy = (

    free_energy_masked[
        np.isfinite(
            free_energy_masked
        )
    ]

)


if len(finite_energy) == 0:

    raise RuntimeError(
        "No valid free-energy values remain "
        "after probability filtering."
    )


# ============================================================
# ACTUAL FREE-ENERGY RANGE
# ============================================================

actual_minimum = np.min(
    finite_energy
)


actual_maximum = np.max(
    finite_energy
)


# ============================================================
# DETERMINE DISPLAY MAXIMUM
# ============================================================

if DISPLAY_MODE == "percentile":

    display_max = np.percentile(

        finite_energy,

        ENERGY_PERCENTILE

    )


elif DISPLAY_MODE == "fixed":

    display_max = MAX_FREE_ENERGY


else:

    raise ValueError(
        "DISPLAY_MODE must be "
        "'percentile' or 'fixed'."
    )


# ============================================================
# SAFETY CHECK
# ============================================================

display_max = max(

    display_max,

    actual_minimum + 0.001

)


# ============================================================
# CLIP ONLY FOR VISUALIZATION
# ============================================================

free_energy_plot = (
    free_energy_masked.copy()
)


free_energy_plot[
    free_energy_plot > display_max
] = display_max


# ============================================================
# PRINT FREE-ENERGY INFORMATION
# ============================================================

print()
print("=" * 70)
print(" FREE-ENERGY ANALYSIS")
print("=" * 70)

print(
    f"Minimum free energy: "
    f"{actual_minimum:.4f} kJ/mol"
)

print(
    f"Maximum calculated free energy "
    f"(sampled region): "
    f"{actual_maximum:.4f} kJ/mol"
)

print(
    f"Display maximum: "
    f"{display_max:.4f} kJ/mol"
)

print(
    f"Probability threshold: "
    f"{MIN_PROBABILITY_FRACTION:.4f} × Pmax"
)

print(
    f"Sampled bins: "
    f"{np.sum(sampled_mask)} / "
    f"{sampled_mask.size}"
)

print(
    f"Smoothing: σ = "
    f"{SMOOTHING:.2f}"
)


# ============================================================
# SAVE NUMERICAL DATA
# ============================================================

"""
CSV columns:

    PC1
    PC2
    Probability
    Free_Energy_kJ_mol
    Sampled

IMPORTANT:

The Free_Energy_kJ_mol column contains the ORIGINAL
calculated values.

It is NOT clipped to display_max.
"""

df = pd.DataFrame({

    "PC1":
        X.flatten(),

    "PC2":
        Y.flatten(),

    "Probability":
        probability.flatten(),

    "Free_Energy_kJ_mol":
        free_energy.flatten(),

    "Sampled":
        sampled_mask.flatten()

})


df.to_csv(

    OUTPUT_CSV,

    index=False

)


print()
print(
    "Numerical data saved:"
)

print(
    "    ",
    OUTPUT_CSV
)


# ============================================================
# CREATE FIGURE
# ============================================================

fig = plt.figure(

    figsize=(7.2, 9.5)

)


# ============================================================
# COLOR MAP
# ============================================================

"""
Use the classic jet-like color scheme seen in the
reference figure:

    low energy:
        dark blue

    intermediate:
        cyan / green / yellow

    high energy:
        orange / red
"""

cmap = plt.colormaps["jet"]


# ============================================================
# COLOR NORMALIZATION
# ============================================================

norm = Normalize(

    vmin=0.0,

    vmax=display_max

)


# ============================================================
# 3D FREE-ENERGY LANDSCAPE
# ============================================================

ax1 = fig.add_subplot(

    2,
    1,
    1,

    projection="3d"

)


# ============================================================
# PREPARE 3D SURFACE
# ============================================================

surface_Z = (
    free_energy_plot.copy()
)


"""
Poorly sampled regions are represented by NaN.

For visualization of the 3D landscape, these regions are
raised to the maximum displayed energy.

This produces the characteristic high-energy plateau
around poorly sampled regions, similar to the reference
shape.
"""

surface_Z[
    ~np.isfinite(surface_Z)
] = display_max


# ============================================================
# 3D SURFACE
# ============================================================

surface = ax1.plot_surface(

    X,
    Y,
    surface_Z,

    cmap=cmap,

    norm=norm,

    linewidth=0,

    antialiased=True,

    rstride=1,

    cstride=1,

    shade=False

)


# ============================================================
# 3D AXIS LABELS
# ============================================================

ax1.set_xlabel(

    "PC1",

    fontsize=10,

    labelpad=2

)


ax1.set_ylabel(

    "PC2",

    fontsize=10,

    labelpad=2

)


ax1.set_zlabel(

    "Free Energy (kJ/mol)",

    fontsize=10,

    labelpad=5

)


# ============================================================
# 3D VIEW
# ============================================================

ax1.view_init(

    elev=ELEVATION,

    azim=AZIMUTH

)


# ============================================================
# 3D Z LIMIT
# ============================================================

ax1.set_zlim(

    0,

    display_max

)


# ============================================================
# 3D TICKS
# ============================================================

ax1.tick_params(

    axis="both",

    labelsize=8

)

ax1.tick_params(

    axis="z",

    labelsize=8

)


# ============================================================
# 3D COLORBAR
# ============================================================

cbar1 = fig.colorbar(

    surface,

    ax=ax1,

    shrink=0.58,

    pad=0.08,

    aspect=20

)


cbar1.set_label(

    "Free Energy (kJ/mol)",

    fontsize=9

)


cbar1.ax.tick_params(

    labelsize=8

)


# ============================================================
# 2D FREE-ENERGY LANDSCAPE
# ============================================================

ax2 = fig.add_subplot(

    2,
    1,
    2

)


# ============================================================
# MASK INVALID VALUES
# ============================================================

plot_2d = np.ma.masked_invalid(

    free_energy_plot

)


# ============================================================
# CONTOUR LEVELS
# ============================================================

levels = np.linspace(

    0,

    display_max,

    N_LEVELS

)


# ============================================================
# FILLED CONTOUR
# ============================================================

filled = ax2.contourf(

    X,
    Y,

    plot_2d,

    levels=levels,

    cmap=cmap,

    vmin=0,

    vmax=display_max,

    extend="max"

)


# ============================================================
# CONTOUR LINES
# ============================================================

"""
The reference uses visible contour lines.

Red contour lines are retained here to make the
free-energy basins easy to identify.
"""

ax2.contour(

    X,
    Y,

    plot_2d,

    levels=levels,

    colors="red",

    linewidths=0.45,

    alpha=0.95

)


# ============================================================
# 2D AXIS LABELS
# ============================================================

ax2.set_xlabel(

    "PC1",

    fontsize=10

)


ax2.set_ylabel(

    "PC2",

    fontsize=10

)


# ============================================================
# 2D TICKS
# ============================================================

ax2.tick_params(

    labelsize=8

)


# ============================================================
# 2D COLORBAR
# ============================================================

cbar2 = fig.colorbar(

    filled,

    ax=ax2,

    pad=0.02,

    aspect=25

)


cbar2.set_label(

    "Free Energy (kJ/mol)",

    fontsize=9

)


cbar2.ax.tick_params(

    labelsize=8

)


# ============================================================
# FIGURE LAYOUT
# ============================================================

plt.subplots_adjust(

    left=0.10,

    right=0.88,

    top=0.97,

    bottom=0.07,

    hspace=0.20

)


# ============================================================
# SAVE PNG
# ============================================================

plt.savefig(

    OUTPUT_PNG,

    dpi=600,

    bbox_inches="tight",

    facecolor="white"

)


# ============================================================
# SAVE PDF
# ============================================================

plt.savefig(

    OUTPUT_PDF,

    dpi=600,

    bbox_inches="tight",

    facecolor="white"

)


# ============================================================
# FINAL INFORMATION
# ============================================================

print()
print("=" * 70)
print(" ANALYSIS COMPLETED SUCCESSFULLY")
print("=" * 70)
print()

print(
    "PNG :",
    OUTPUT_PNG
)

print(
    "PDF :",
    OUTPUT_PDF
)

print(
    "CSV :",
    OUTPUT_CSV
)

print()

print(
    "Free-energy equation:"
)

print(
    "    ΔG = -RT ln(P/Pmax)"
)

print(
    f"Temperature: "
    f"{TEMPERATURE:.1f} K"
)

print(
    f"Histogram: "
    f"{NBINS} × {NBINS}"
)

print(
    f"Gaussian smoothing: "
    f"σ = {SMOOTHING:.2f}"
)

print(
    f"Display mode: "
    f"{DISPLAY_MODE}"
)

print(
    f"Display maximum: "
    f"{display_max:.4f} kJ/mol"
)

print()

print(
    "Visualization:"
)

print(
    "    Colormap: jet"
)

print(
    "    Low energy: blue"
)

print(
    "    High energy: red"
)

print()

print(
    "Closing figure..."
)


# ============================================================
# DISPLAY
# ============================================================

plt.show()