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

The main code (`Cosmic String Extraction Statistics.py`) constructs a simulated environment for detecting cosmic string wakes.

### 1. Initialize Universe and Parameters:
Sets Hubble volume, cosmological parameters, and random seeds for reproducibility.

```python
hubbleParameter = 70.41 * 1000   # [m/s/Mpc] Hubble constant H(z=0)
densityRatio = 0.3               # Non-relativistic matter density for flat ΛCDM
speedOfLight = 299792458         # [m/s] Speed of light
hubbleScale = speedOfLight / hubbleParameter  # [Mpc] Hubble length
```

### 2. Cosmic String Wake Generation:
Builds a finite-length cosmic string segment, determines its wake geometry, applies rotations from $SO(3)$ group elements to orient it arbitrarily in space, and then maps it into redshift space using distance-redshift relations.

```python
Gμ = 3E-7  # string tension
deficitAngle = 8 * np.pi * Gμ
wakeLength = c1 * formationTime * speedOfLight * meterToMpc  # [Mpc]
wakeDepth = formationTime * gammaFactor * stringSpeed * meterToMpc  # [Mpc] radial length
```

The wake geometry is defined by six vertices:

```python
wakeWedge = np.array([
    wakeTipPoint, wakeEndPoints[0], wakeEndPoints[1],
    projectedWakePoints[0], projectedWakePoints[1], projectedWakePoints[2]
])
```

### 3. Wake in Physical Space
Assigning Temperature Fields:
Within the wake's convex hull, assigns a temperature gradient based on the brightness temperature formula $\delta T_b(\nu)$. Outside the wake, the temperature field is ambient.

```python
def brightnessTemperature(z):
    deexcitationCrossSection = 0.16
    if (kineticTemperature(z) > 3*gasTemperature(z)):
        u = kineticTemperature(z)
    else:
        u = 3*gasTemperature(z)
    if (kineticTemperature(z) > 3*gasTemperature(z)):
        c = 4
    else:
        c = 1 + kineticTemperature(z)/gasTemperature(z)
    a = 0.017*deexcitationCrossSection/(1+deexcitationCrossSection)*(1-photonCMBTemp(z)/u)*np.sqrt(1+z)*c
    return a
```

### 4. Primordial Noise Embedding:
Loads or simulates a 3D $\Lambda$CDM cosmological noise map (using 21cmFAST simulations) and superimposes it with the wake signal. This results in a realistic data cube containing both signal and noise.

```python
perturbationNoise = box * 0.001  # converting from mK to K
smallNoiseMap = np.copy(thirdDimensionSlice)  # 29x29x29 array of temps
reshapedTempMap = gridTemps.reshape((29,29,29))  # 29x29x29 array of temps
combinedMap = smallNoiseMap + reshapedTempMap
```

### 5. Matched Filtering & Statistics:
Performs matched filtering by correlating the wake template with the data. Evaluates various slices and orientations, unfolding 2D and 3D data arrays into 1D arrays if necessary to maximize the signal-to-noise ratio.

```python
def unfolder(grid, type):
    amplitudeValues = []
    virtualGrid = np.copy(grid)
    if type == 'horizontal':
        for row in virtualGrid:
            for value in row:
                amplitudeValues.append(value)
    elif type == 'vertical':
        for column in virtualGrid.T:
            for value in column:
                amplitudeValues.append(value)
    return amplitudeValues

horizontalUnfoldedWake = unfolder(convolvedWake, 'horizontal')
verticalUnfoldedWake = unfolder(convolvedWake, 'vertical')

The matched filtering is then applied:

python
horizontalWakeWake = np.correlate(horizontalUnfoldedConvolvedWake, horizontalUnfoldedConvolvedWake, mode='full')
horizontalWakeNoise = np.correlate(horizontalUnfoldedConvolvedWake, horizontalUnfoldedConvolvedNoise, mode='full')
horizontalWakeCombined = np.correlate(horizontalUnfoldedConvolvedWake, horizontalUnfoldedConvolvedCombined, mode='full')
```



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


## Theoretical Background

When a scalar field $\phi$ with a potential

\[
V(\phi) = \frac{1}{4}g(|\phi|^2 - \eta^2)^2
\]

undergoes spontaneous symmetry breaking, the vacuum manifold $M$ can become topologically non-trivial. This occurs when the scalar field transitions to a degenerate vacuum state below a critical temperature $T_c \sim \eta$, forming a manifold of vacua $M \simeq S^1$. 

**Key Insight:** The first homotopy group of the vacuum manifold, $\pi_1(M)$, is non-trivial ($\pi_1(S^1) = \mathbb{Z}$), ensuring the existence of stable one-dimensional line defects — cosmic strings. These defects form due to discontinuities in the field configuration along closed loops, which cannot be continuously shrunk to a point.

The resulting cosmic strings generate planar overdensities (wakes) in the surrounding matter. These overdensities perturb the cosmic microwave background (CMB) and induce temperature anisotropies observable in 21cm maps. The brightness temperature fluctuation $\delta T_b(\nu)$ at a redshift $z$ is given by:

\[
\delta T_b (\nu) = [0.07 \, \text{K}] \frac{x_c}{1+x_c} \left(1 - \frac{T_\gamma(z)}{T_{K/g}(z)}\right) \sqrt{1+z},
\]

where:

- $x_c$: Collision coupling coefficient.
- $T_\gamma(z)$: Temperature of CMB photons.
- $T_{K/g}(z)$: Kinetic (gas) temperature of the wake.

**Visualization of Symmetry Breaking Potential:**  
\includegraphics[scale=0.5]{https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/GoldstonePotential.png}

Above: The scalar field’s vacuum manifold $S^1$ emerges as the minimum of the symmetry-breaking potential $V(\phi)$, creating a topologically non-trivial space.

---

### Cosmic String Wakes in Redshift Space

Cosmic strings moving relativistically at velocity $v_s$ induce overdensities (wakes) with a conical geometry in physical space. The wake’s geometry is defined by the deficit angle:

\[
\alpha = 8 \pi G \mu,
\]

where:

- $G\mu$: Dimensionless string tension.
- $\mu$: Energy per unit length of the string.

To represent the wake in redshift space, the cosmological distance-redshift relation is used:

\[
d(z) = \int_0^z \frac{c \, dz'}{H(z')},
\]

where:

- $c$: Speed of light.
- $H(z)$: Hubble parameter.

The wake is narrow in redshift space due to the small value of the deficit angle and the localized nature of its temperature gradient $\delta T_b(\nu)$.

**Figure: Wake Geometry in Redshift Space**  
\includegraphics[scale=0.5]{https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/wakeRedshiftSpace.png}

Above: A cosmic string wake mapped into redshift space, showing its finite extent in the transverse and longitudinal directions. The wake's narrow geometry in the redshift axis is key to its detection.

---

### Matched Filtering: Extracting the Signal

The brightness temperature fluctuation $\delta T_b (\nu)$ is buried under primordial $\Lambda$CDM perturbations, making detection challenging. Matched filtering is employed to extract the signal by correlating the observed data with a theoretical template (the wake profile). 

The matched filtering statistic is given by:

\[
s(t) = \sum_k h(t-k) d(k),
\]

where:

- $s(t)$: Match filter amplitude.
- $h(k)$: Wake template.
- $d(k)$: Observed data.

This convolution maximizes the signal-to-noise ratio, making it possible to identify subtle wake signals in noisy maps.

**Figure: 1D Matched Filtering**  
\includegraphics[scale=0.5]{https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/matchFilter.png}

Above: Match filter results for various data sets. The peak indicates the presence of a cosmic string wake in noisy data, confirming its detection.




> [!TIP]
> For a detailed derivation of the brightness temperature formula, the homotopy classification of defects, and scaling solutions for cosmic strings, consult the PDF in the main repository.

>[!NOTE]
> Once the method is validated on simulated data, it can be applied to actual 21cm observations (e.g., from upcoming radio interferometers) to search for cosmic string imprints in the real universe.
