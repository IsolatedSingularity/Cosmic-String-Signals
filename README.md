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

<details>
  <summary><i>Show Detailed Python Code from the Provided Repository</i></summary>

```python
#%% Importing modules 유
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from astropy.cosmology import Planck18
from astropy.cosmology import z_at_value
import astropy.units as u
import random as rd
from scipy import signal

#%% Defining compact plotting functions 유

def quickPlot(scale,points,labels,colors,colormaps,convexhull,testPoints=[]):
    fig = plt.figure(figsize=(12, 10), dpi=80)
    ax = fig.add_subplot(111, projection="3d")
    approximateCenter = []
    for i in points:
        if len(i) == 6:
            for j in range(0,3):
                approximateCenter.append(
                    0.5*(np.max(i.T[j])-np.min(i.T[j])) + np.min(i.T[j])
                )
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('y [Mpc]')
    ax.set_zlabel('z [Mpc]')
    ax.set_xlim(approximateCenter[0]-scale/2.7,approximateCenter[0]+scale/2.7)
    ax.set_ylim(approximateCenter[1]-scale/2.7,approximateCenter[1]+scale/2.7)
    ax.set_zlim(approximateCenter[2]-scale/2.7,approximateCenter[2]+scale/2.7)
    
    if convexhull == False:
        if len(colormaps) != 0:
            physicalPlot = ax.scatter(
                points[-1][:,0],points[-1][:,1],points[-1][:,2], c = colormaps,
                cmap = plt.cm.viridis, s = 90
            )
            fig.colorbar(physicalPlot, ax=ax, label = 'δT(z)', pad=0.1)
        if len(colormaps) == 0:
            for pointType in range(len(points)):
                if points[pointType].ndim == 2:
                    ax.scatter(
                        points[pointType][:,0], points[pointType][:,1], points[pointType][:,2],
                        label = labels[pointType], color = colors[pointType], s = 90
                    )
                else:
                    ax.scatter(
                        points[pointType][0], points[pointType][1], points[pointType][2],
                        label = labels[pointType], color = colors[pointType], s = 90
                    )
    if convexhull != False:
        hull = ConvexHull(points[0])
        ax.scatter(
            points[0].T[0], points[0].T[1], points[0].T[2], color='red', s = 90
            )
        for s in hull.simplices:
            s = np.append(s, s[0]) #simplice
            ax.plot(points[0][s, 0], points[0][s, 1], points[0][s, 2], color='maroon')
        if len(testPoints) != 0:
            ax.scatter(
            testPoints[:,0],testPoints[:,1],testPoints[:,2], color='black', s = 90
            )
    ax.legend()
    plt.tight_layout()
    plt.autoscale()
    plt.show()
    return
```
# Additional plotting, rotation, scaling, wedge construction, redshift conversion,
# matched filtering, and convolution functions are defined below. 
# For brevity, we omit re-displaying all of them here, but they appear in the full code.

</details> <details> <summary><i>Show LaTeX Theory Snippet from the Report</i></summary>

\section{Cosmic String Wakes and Signal Extraction}

Topological defects arise when the vacuum manifold $\mathcal{M}$ of a field theory is non-trivial.
For a $U(1)$ symmetry breaking, the vacuum manifold after symmetry breaking is isomorphic to
$M \simeq U(1) \simeq S^1$. The non-trivial fundamental group $\pi_1(S^1)=\mathbb{Z}$ guarantees
the existence of stable line defects (cosmic strings).

The tension of such a string is given by:
\begin{equation}
G\mu \sim \left(\frac{\eta}{M_{\mathrm{Pl}}}\right)^2,
\end{equation}
and the resulting wake in a 21cm map produces a brightness temperature shift:
\begin{equation}
\delta T_b(\nu) = [0.07 \, K]\frac{x_c}{1+x_c}\left(1-\frac{T_\gamma(z)}{T_{K/g}(z)}\right)\sqrt{1+z}.
\end{equation}

To convert physical distances to redshifts:
\begin{equation}
d(z)=\int_0^z \frac{c}{H(z')}dz',
\end{equation}
and we use numerical routines to invert this relation, embedding the cosmic string wake into
redshift space.

Matched filtering $s(t)$ for a 1D unfolded data set is given by:
\begin{equation}
s(t) = \sum_{k=-n^2/2}^{n^2/2} h[t-k]d[k],
\end{equation}
where $h$ is the template (e.g., our modeled wake) and $d$ is the data (wake+noise).

</details>
Caveats

    Convex Hull Limitations:
    When the actual deficit angle $\alpha = 8\pi G\mu$ is tiny, the wake is almost a plane, making convex hull detection challenging. The current algorithm relies on simplices, which fails in near-1D geometries.

    Redshift Conversion Issues:
    Using astropy.cosmology functions, extremely small or large $z$ values may cause numerical convergence problems. Thus, the code currently focuses on realistic $z$ ranges and approximates comoving coordinates.

    Simplifications in Physics:
    The temperature model $\delta T_b(\nu)$ assumes certain simplifications about gas thermodynamics and the kinetic temperature relationship. Real-world complexities (shock heating, non-linear structure) may require more sophisticated modeling.

Next Steps

Implement match filtering on multiple slicing orientations to find the most robust direction for signal extraction.
Explore alternative geometric detection algorithms that handle planar or line-like topologies robustly.
Use MCMC or Bayesian inference frameworks with the matched filtering outputs to place statistical constraints on $G\mu$.

    Integrate machine learning or wavelet-based analysis for non-Gaussian features in the noise field.

    [!TIP] For a detailed derivation of the brightness temperature formula, the homotopy classification of defects, and scaling solutions for cosmic strings, consult the PDF in the main repository.

    [!NOTE] Once the method is validated on simulated data, it can be applied to actual 21cm observations (e.g., from upcoming radio interferometers) to search for cosmic string imprints in the real universe.
