# Topological-Defect-Signals
###### Under the supervision of [Professor Robert Brandenberger](https://www.physics.mcgill.ca/~rhb/) at [McGill University](https://www.mcgill.ca/) and in collaboration with [Mattéo Blamart](https://inspirehep.net/authors/2077637). Funded by the [NSERC USRA](https://www.nserc-crsng.gc.ca/students-etudiants/ug-pc/usra-brpc_eng.asp) award with [FRQNT supplements](https://frq.gouv.qc.ca/en/program/supplements-of-the-nserc-undergraduate-student-research-awards-usra-bpca-2023-2024/).

![alt text](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/2DConvolution.png?raw=true)

## Objective

In the very early universe, symmetry-breaking phase transitions can give rise to non-trivial vacuum manifolds. These manifolds host **topological defects**, whose stability and existence are governed by the homotopy groups of the vacuum manifold. Among these defects, **cosmic strings** are one-dimensional line defects associated with a spontaneously broken U(1) symmetry.

The energy scale at which the U(1) symmetry is broken sets the tension $G\mu$ of the cosmic strings, where $G$ is Newton’s gravitational constant and $\mu$ is the mass per unit length of the string. The tension is related to the symmetry-breaking scale $\eta$ by an approximate quadratic relation:

$$
G\mu \sim \left(\frac{\eta}{M_{\text{Pl}}}\right)^2,
$$

where $M_{\text{Pl}}$ is the Planck mass. After their formation, cosmic strings stretch across the universe and as they move at relativistic speeds, they generate **wakes** in the surrounding matter distribution. These wakes induce temperature anisotropies observable in 21cm radiation. The challenge is that these signals are deeply buried under primordial $\Lambda$CDM perturbations, which act as a non-linear form of noise.

**Goal:** Develop advanced statistical methods, such as matched filtering and correlation analyses, to extract the cosmic string wake signal from this noisy background. By doing so, we can potentially identify or constrain the presence of cosmic strings, providing insights into high-energy physics beyond the Standard Model and shedding light on the early universe’s structure formation.

## Theoretical Background

When a scalar field $\phi$ with a potential

$$
V(\phi) = \frac{1}{4}g(|\phi|^2 - \eta^2)^2
$$

undergoes spontaneous symmetry breaking, the vacuum manifold $M$ can become topologically non-trivial. If $M \simeq S^1$, the first homotopy group $\pi_1(M)$ is non-trivial, ensuring the existence of stable line defects — cosmic strings.

These strings create planar overdensities (wakes) as they pass through matter. The resulting brightness temperature fluctuation at redshift $z$ in 21cm maps can be modeled as:

$$
\delta T_b(\nu) = [0.07 \, \text{K}] \, \frac{x_c}{1+x_c} \left(1 - \frac{T_{\gamma}(z)}{T_{K/g}(z)}\right) \sqrt{1+z},
$$

where $x_c$ is the collision coupling coefficient, $T_{\gamma}(z)$ is the CMB photon temperature, and $T_{K/g}(z)$ characterizes the kinetic/gas temperature within the wake.

To convert physical coordinates into redshift space, we solve the cosmological distance-redshift relation. For a flat FRW universe:

$$
d(z) = \int_0^z \frac{c \, dz'}{H(z')},
$$

where $H(z)$ is the Hubble parameter. Numerical methods (e.g., `astropy.cosmology`) provide this mapping, allowing the wake—originally defined in physical coordinates—to be represented in redshift space, where the anisotropies are more directly comparable to observations.

However, small angles and subtle temperature differences demand sophisticated statistics to extract the signal from noise. **Matched filtering**, a technique widely used in signal processing (e.g., gravitational wave detections), is employed here. By correlating an assumed template (wake profile) with the observed data, we amplify the signal relative to the random noise.

## Code Functionality

The main code (`Cosmic String Extraction Statistics.py`) constructs a simulated environment:

1. **Initialize Universe and Parameters:**  
   Sets Hubble volume, cosmological parameters, and random seeds for reproducibility.

2. **Cosmic String Wake Generation:**  
   Builds a finite-length cosmic string segment, determines its wake geometry, applies rotations from $SO(3)$ group elements to orient it arbitrarily in space, and then maps it into redshift space using distance-redshift relations.

3. **Assigning Temperature Fields:**  
   Within the wake’s convex hull, assigns a temperature gradient based on the brightness temperature formula $\delta T_b(\nu)$. Outside the wake, the temperature field is ambient.

4. **Primordial Noise Embedding:**  
   Loads or simulates a 3D $\Lambda$CDM cosmological noise map (using 21cmFAST simulations) and superimposes it with the wake signal. This results in a realistic data cube containing both signal and noise.

5. **Matched Filtering & Statistics:**  
   Performs matched filtering by correlating the wake template with the data. Evaluates various slices and orientations, unfolding 2D and 3D data arrays into 1D arrays if necessary to maximize the signal-to-noise ratio. Also explores wavelet transforms and other correlation methods to further isolate the signal.


</details>

## Caveats

- Convex Hull Limitations:
When the actual deficit angle $\alpha = 8\pi G\mu$ is tiny, the wake is almost a plane, making convex hull detection challenging. The current algorithm relies on simplices, which fails in near-1D geometries.

- Redshift Conversion Issues:
Using astropy.cosmology functions, extremely small or large $z$ values may cause numerical convergence problems. Thus, the code currently focuses on realistic $z$ ranges and approximates comoving coordinates.

- Simplifications in Physics:
The temperature model $\delta T_b(\nu)$ assumes certain simplifications about gas thermodynamics and the kinetic temperature relationship. Real-world complexities (shock heating, non-linear structure) may require more sophisticated modeling.

## Next Steps

- [x] Implement match filtering on multiple slicing orientations to find the most robust direction for signal extraction.
- [ ] Explore alternative geometric detection algorithms that handle planar or line-like topologies robustly.
- [ ] Use MCMC or Bayesian inference frameworks with the matched filtering outputs to place statistical constraints on $G\mu$.
- [ ] Integrate machine learning or wavelet-based analysis for non-Gaussian features in the noise field.

> [!TIP]
> For a detailed derivation of the brightness temperature formula, the homotopy classification of defects, and scaling solutions for cosmic strings, consult the PDF in the main repository.

>[!NOTE]
> Once the method is validated on simulated data, it can be applied to actual 21cm observations (e.g., from upcoming radio interferometers) to search for cosmic string imprints in the real universe.


Under the supervision of Professor Robert Brandenberger at McGill University, and in collaboration with Mattéo Blamart. Funded by the NSERC USRA award with FRQNT supplements.
alt text
Objective
Phase transitions in the early universe caused spontaneous symmetry breaking events, leading to the formation of topological defects. Linear topological defects associated with the U(1) symmetry group are known as cosmic strings, which create high-energy signals in the form of wakes as they propagate through spacetime. These wakes leave observable imprints, such as temperature anisotropies in the cosmic microwave background (CMB), making them critical for studying high-energy physics and early universe structure formation. The purpose of this repository is to:

    Model finite-length cosmic string wakes in physical and redshift space.
    Embed these wakes in noisy cosmological maps.
    Extract the wake signal using match filtering and statistical techniques.

Code Functionality
Wake Construction in Physical Space
Cosmic strings create planar overdensities or "wakes" as they propagate through spacetime. The geometry of these wakes is determined by the deficit angle:
α=8πGμ,
α=8πGμ,
where GμGμ is the dimensionless string tension. The wake is modeled as a wedge-shaped structure defined by six vertices. Below is an example Python snippet illustrating this construction:

python
import numpy as np

# Define wake parameters
G_mu = 3e-7  # String tension
deficit_angle = 8 * np.pi * G_mu
wake_length = 100  # Mpc
wake_depth = 10     # Mpc

# Define vertices of the wedge
vertices = np.array([
    [0, 0, 0],
    [wake_depth * np.cos(deficit_angle / 2), wake_depth * np.sin(deficit_angle / 2), 0],
    [wake_depth * np.cos(deficit_angle / 2), -wake_depth * np.sin(deficit_angle / 2), 0],
    [0, 0, wake_length],
    [wake_depth * np.cos(deficit_angle / 2), wake_depth * np.sin(deficit_angle / 2), wake_length],
    [wake_depth * np.cos(deficit_angle / 2), -wake_depth * np.sin(deficit_angle / 2), wake_length]
])

Visualization:
The wake geometry in physical space is visualized below: alt text
Redshift Conversion and Scaling
To match cosmological observations, we rescale the vertical axis (in Mpc) to redshift space using a numerical relationship between comoving distance and redshift derived from FRW cosmology:

python
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u

# Convert distance (Mpc) to redshift
distance_to_redshift = lambda d: cosmo.z_at_value(cosmo.comoving_distance, d * u.Mpc)
redshift_to_distance = lambda z: cosmo.comoving_distance(z).value

# Example usage
z = distance_to_redshift(100)  # Redshift corresponding to 100 Mpc
d = redshift_to_distance(1)    # Distance corresponding to z=1

Visualization:
The scaling relationship between distance and redshift is shown below: alt text
Signal Extraction via Match Filtering
Match filtering is employed to extract the cosmic string signal from noisy maps. For a one-dimensional dataset, it is defined as:
s(t)=∑kh[t−k]⋅d[k],
s(t)=k∑​h[t−k]⋅d[k],
where:

    s(t)s(t): Match filter amplitude,
    h(t)h(t): Template signal (e.g., cosmic string wake),
    d(t)d(t): Observed data (signal + noise).

Implementation Example:

python
from scipy.signal import correlate

# Define template (wake signal) and data (signal + noise)
template = np.array([...])   # Wake gradient map
data = np.array([...])       # Combined map with noise

# Compute match filter amplitude
match_filter_amplitude = correlate(data, template, mode='full')

Visualization:
Below are slices of combined maps showing the embedded wake signal: alt text
Results and Figures
Wake Geometry in Redshift Space
After rescaling and applying transformations for an expanding universe, the wake geometry appears as follows: alt text
Match Filtering Results
Match filtering results for different datasets demonstrate how the cosmic string signal can be extracted from noisy maps: alt text
Caveats

    Convex Hull Algorithm: The convex hull approach struggles when the deficit angle αα is small (α=8πGμα=8πGμ), causing the wake to approximate a plane.
    Scaling Issues: The conversion between physical and comoving coordinates relies on numerical approximations that may not converge at extreme values of redshift.

Next Steps
At this stage, we've clarified the problem statement and established a foundational code framework. Future work includes:

    Generalizing algorithms for higher-dimensional topological defect signals.
    Developing an algorithm to invert the redshift function without convergence limitations.
    Enhancing signal extraction methods using machine learning techniques.

This README now incorporates detailed descriptions, equations, figures, and code snippets while maintaining a structure similar to your desired format.
