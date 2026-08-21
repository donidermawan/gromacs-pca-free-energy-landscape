# GROMACS-PCA-free-energy-landscape

A Python script for generating **2D and 3D Principal Component Analysis (PCA) free-energy landscapes** from GROMACS `.xvg` PCA data.

The script converts PC1/PC2 coordinates into a 2D probability distribution, optionally applies light Gaussian smoothing, and calculates the relative free energy using:

```text
ΔG = -RT ln(P / Pmax)
```

where the most populated region is defined as approximately **0 kJ/mol**.

The resulting landscape uses a **blue → cyan → green → yellow → orange → red** color scheme, with low-energy regions represented by blue and high-energy regions represented by red.

---

## Features

* Reads GROMACS `.xvg` files containing PC1 and PC2 coordinates
* Generates a 2D PC1/PC2 population histogram
* Normalizes the probability distribution
* Calculates relative free energy in kJ/mol
* Produces a **3D free-energy surface**
* Produces a **2D filled free-energy contour map**
* Identifies poorly sampled regions using a probability threshold
* Allows automatic or fixed free-energy display limits
* Exports high-resolution PNG figures at 600 DPI
* Exports PDF figures suitable for publications and reports
* Exports the numerical PC1/PC2 probability and free-energy grid as CSV

---

## Requirements

Python 3 with the following packages:

```bash
pip install numpy pandas matplotlib scipy
```

The required packages are:

* `numpy`
* `pandas`
* `matplotlib`
* `scipy`

---

## Input

The default input file is:

```text
pc1_pc2_2d.xvg
```

The file should contain at least two numerical columns:

```text
PC1    PC2
```

For example:

```text
-2.1345    1.2456
-2.1021    1.2874
-2.0876    1.3012
-2.0543    1.3421
...
```

GROMACS metadata lines beginning with:

```text
@
#
```

are automatically ignored.

Inline comments beginning with `#` are also ignored.

---

## Usage

Place the Python script and the input file in the same directory:

```text
PCA_Free_Energy_Landscape.py
pc1_pc2_2d.xvg
```

Run:

```bash
python3 PCA_Free_Energy_Landscape.py
```

The script will read the PCA trajectory, calculate the probability distribution and relative free-energy surface, generate the 3D and 2D landscapes, and save the numerical results.

---

## Main Parameters

The main analysis parameters can be modified directly in the script:

```python
INPUT_FILE = "pc1_pc2_2d.xvg"

TEMPERATURE = 300.0

NBINS = 100

SMOOTHING = 0.5

N_LEVELS = 20

DISPLAY_MODE = "percentile"

ENERGY_PERCENTILE = 95

MAX_FREE_ENERGY = 10.0

MIN_PROBABILITY_FRACTION = 0.001
```

---

## Parameters

| Parameter                  | Description                                                                    |          Default |
| -------------------------- | ------------------------------------------------------------------------------ | ---------------: |
| `INPUT_FILE`               | GROMACS PCA `.xvg` input file                                                  | `pc1_pc2_2d.xvg` |
| `TEMPERATURE`              | Simulation temperature in K                                                    |          `300.0` |
| `NBINS`                    | Number of histogram bins in each PCA dimension                                 |            `100` |
| `SMOOTHING`                | Gaussian smoothing parameter                                                   |            `0.5` |
| `N_LEVELS`                 | Number of contour levels                                                       |             `20` |
| `DISPLAY_MODE`             | Method used to determine the displayed maximum free energy                     |     `percentile` |
| `ENERGY_PERCENTILE`        | Percentile used for the displayed maximum when using percentile mode           |             `95` |
| `MAX_FREE_ENERGY`          | Fixed displayed maximum when using fixed mode                                  |           `10.0` |
| `MIN_PROBABILITY_FRACTION` | Minimum probability relative to `Pmax` used to identify poorly sampled regions |          `0.001` |
| `ELEVATION`                | 3D viewing elevation angle                                                     |             `28` |
| `AZIMUTH`                  | 3D viewing azimuth angle                                                       |           `-125` |

---

## Histogram Resolution

The PC1/PC2 trajectory is converted into a 2D histogram.

The default is:

```python
NBINS = 100
```

which produces:

```text
100 × 100
```

histogram bins.

Increasing `NBINS` provides higher spatial resolution but can increase noise when the trajectory contains insufficient sampling.

Decreasing `NBINS` produces a smoother and more generalized landscape.

For most typical PCA trajectories, `100 × 100` is a useful starting point.

---

## Gaussian Smoothing

The script applies optional Gaussian smoothing to the histogram before calculating the free energy.

The default is:

```python
SMOOTHING = 0.5
```

The general behavior is:

```text
SMOOTHING = 0.0
    No smoothing

SMOOTHING = 0.5
    Very light smoothing

SMOOTHING = 1.0
    Light smoothing

SMOOTHING = 2.0
    Moderate smoothing

SMOOTHING = 3.0
    Strong smoothing
```

A relatively small value is intentionally used because the desired landscape should retain **distinct and relatively sharp free-energy basins** rather than becoming excessively smooth.

If a very sharp, unsmoothed landscape is required, use:

```python
SMOOTHING = 0.0
```

However, completely removing smoothing can increase noise in sparsely sampled regions.

---

## Free Energy Calculation

The probability distribution is normalized and converted into relative free energy according to:

```text
ΔG = -RT ln(P / Pmax)
```

where:

```text
R    = 0.008314462618 kJ mol⁻¹ K⁻¹
T    = simulation temperature in Kelvin
P    = probability of a PC1/PC2 bin
Pmax = maximum probability
```

Because the probability is normalized relative to the maximum probability:

```text
P = Pmax
```

gives:

```text
ΔG = 0 kJ/mol
```

Therefore, the most populated region of the PCA space corresponds to the minimum relative free energy.

The calculated values are **relative free energies**, not absolute thermodynamic free energies.

---

## Probability Threshold

Very poorly populated histogram bins can produce extremely large or undefined free-energy values.

The script therefore uses:

```python
MIN_PROBABILITY_FRACTION = 0.001
```

This means that bins with:

```text
P < 0.001 × Pmax
```

are considered poorly sampled for visualization.

These regions are masked in the 2D landscape and treated as high-energy boundaries in the 3D visualization.

---

## Display Free-Energy Range

The script separates the **calculated free-energy values** from the **visualization range**.

By default:

```python
DISPLAY_MODE = "percentile"
```

and:

```python
ENERGY_PERCENTILE = 95
```

Therefore, the displayed maximum is determined from the 95th percentile of the finite, sampled free-energy values.

Values above this display limit are clipped **only for visualization**.

This prevents a small number of extremely high-energy regions from dominating the color scale and making the low-energy basins difficult to distinguish.

### Percentile mode

```python
DISPLAY_MODE = "percentile"
ENERGY_PERCENTILE = 95
```

The display maximum is automatically determined from the data.

### Fixed mode

Alternatively:

```python
DISPLAY_MODE = "fixed"
MAX_FREE_ENERGY = 10.0
```

will display the landscape using a fixed maximum of:

```text
10.0 kJ/mol
```

---

## Output

The script generates three files.

### 1. PCA Free-Energy Landscape

```text
PCA_Free_Energy_Landscape_Refined.png
```

The PNG contains:

* 3D PCA free-energy surface
* 2D PCA free-energy contour map
* Free-energy color scales
* PC1 and PC2 axes
* Free-energy axis in kJ/mol

The figure is exported at:

```text
600 DPI
```

and is suitable for high-resolution figures, presentations, and manuscripts.

---

### 2. Publication-Quality PDF

```text
PCA_Free_Energy_Landscape_Refined.pdf
```

The PDF contains the same 3D and 2D free-energy landscapes and is suitable for:

* Scientific publications
* Manuscripts
* Reports
* Presentations
* Archival figures

---

### 3. Numerical Data

```text
PCA_Free_Energy_Data_Refined.csv
```

The CSV contains:

```text
PC1
PC2
Probability
Free_Energy_kJ_mol
Sampled
```

For example:

```text
PC1,PC2,Probability,Free_Energy_kJ_mol,Sampled
-3.95,-3.90,0.0000012,18.52,False
-3.95,-3.82,0.0000045,15.08,True
...
```

The CSV can be further processed using:

* Python
* R
* MATLAB
* Origin
* Excel
* Other scientific visualization software

### Important

The `Free_Energy_kJ_mol` column contains the **original calculated free-energy values**.

The values are **not clipped to the visualization maximum**.

---

## Workflow

```text
GROMACS PCA Analysis
        │
        ▼
pc1_pc2_2d.xvg
        │
        ▼
Read PC1 / PC2 Coordinates
        │
        ▼
Remove Invalid Values
        │
        ▼
2D Histogram
        │
        ▼
Gaussian Smoothing
        │
        ▼
Probability Distribution
        │
        ▼
Normalize Probability
        │
        ▼
Free Energy Calculation
        │
        ▼
Probability Threshold
        │
        ├───────────────────┐
        ▼                   ▼
  3D Free-Energy       2D Free-Energy
      Surface             Contour Map
        │                   │
        └─────────┬─────────┘
                  ▼
           PNG / PDF / CSV
```

---

## Interpretation

The resulting PCA free-energy landscape represents the conformational free-energy surface sampled by the molecular dynamics trajectory in PC1–PC2 space.

### Low Free-Energy Regions

Low-energy regions correspond to highly populated regions of the PCA trajectory.

These regions may represent relatively stable or frequently visited conformational states.

### Free-Energy Minima

Distinct minima may indicate different conformational basins or metastable states.

The relative depth of a basin reflects its population within the sampled trajectory.

### High Free-Energy Regions

High-energy regions correspond to less frequently populated areas of the PC1–PC2 conformational space.

These regions may occur between major conformational basins and can represent unfavorable or rarely visited conformations.

### Transitions

Connections between neighboring minima can provide qualitative information about conformational transitions along the dominant principal components.

However, a PCA free-energy landscape alone should not be interpreted as a direct kinetic pathway or transition-state calculation.

For kinetic interpretation, additional analyses such as transition-path analysis, Markov state modeling, or appropriate time-dependent analyses may be required.

---

## Relative Free Energy

The free-energy values generated by this script are **relative free energies**.

The reference point is:

```text
Most populated region = 0 kJ/mol
```

Consequently, the numerical values should be interpreted as differences in free energy relative to the most populated region of the sampled PCA space.

The absolute magnitude and quality of the landscape depend on:

* Trajectory length
* Sampling quality
* Histogram resolution
* Gaussian smoothing
* PCA calculation
* Temperature
* Probability threshold
* Conformational diversity represented in the trajectory

---

## Visualization vs. Numerical Data

This script intentionally separates numerical analysis from visualization.

The workflow is:

```text
Original probability
        │
        ▼
Original free energy
        │
        ├───────────────► CSV
        │                 original values
        │
        ▼
Probability filtering
        │
        ▼
Visualization clipping
        │
        ├───────────────► 3D landscape
        │
        └───────────────► 2D contour map
```

Therefore, the displayed landscape should not be assumed to contain the complete numerical free-energy range.

The CSV should be used whenever the original calculated values are required for quantitative analysis.

---

## Example Configuration

A typical configuration is:

```python
INPUT_FILE = "pc1_pc2_2d.xvg"

TEMPERATURE = 300.0

NBINS = 100

SMOOTHING = 0.5

N_LEVELS = 20

DISPLAY_MODE = "percentile"

ENERGY_PERCENTILE = 95

MAX_FREE_ENERGY = 10.0

MIN_PROBABILITY_FRACTION = 0.001
```

This configuration is designed to produce a relatively sharp free-energy landscape while reducing the visual influence of extremely poorly sampled regions.

---

## Recommended Directory Structure

A simple project directory can be organized as:

```text
gromacs-pca-free-energy-landscape/
│
├── PCA_Free_Energy_Landscape.py
├── pc1_pc2_2d.xvg
├── README.md
│
├── PCA_Free_Energy_Landscape_Refined.png
├── PCA_Free_Energy_Landscape_Refined.pdf
└── PCA_Free_Energy_Data_Refined.csv
```

---

## Citation and Acknowledgement

If this script is used in a scientific publication, presentation, thesis, or report, please cite or acknowledge this repository.

The script is intended as a practical tool for visualization and analysis of PCA-derived conformational free-energy landscapes from GROMACS molecular dynamics simulations.

---

## License

This project is provided for scientific and research use.

Please cite or acknowledge the repository if the script contributes to published research.
