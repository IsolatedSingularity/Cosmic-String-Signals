# Cosmic String Wake Signals

#### Under the supervision of [Professor Robert Brandenberger](https://www.physics.mcgill.ca/~rhb/) at [McGill University](https://www.mcgill.ca/), in collaboration with [Mattéo Blamart](https://inspirehep.net/authors/2077637).  

---

## Objective

Cosmic strings—linear topological defects formed during early universe phase transitions—leave detectable signals in spacetime by creating regions of overdensities known as **wakes**. These cosmic string wakes distort the brightness temperature of surrounding hydrogen gas clouds, which can be detected through **21 cm cosmology**.  

However, detecting these signals is challenging due to primordial **ΛCDM noise**. This repository models cosmic string wakes and develops statistical methods, such as **matched filtering**, to extract the wake's faint signal amidst the strong noise fluctuations.

---

## Code Functionality

### **1. Constructing Wakes in 3D Space**

#### Theoretical Background
A **cosmic string wake** is formed as the string propagates through plasma, producing a thin triangular planar overdensity. The wake's defining geometry is rooted in the **deficit angle** caused by the string's gravitational lensing:

$$
\alpha = 8\pi G\mu,
$$


where:
- $$G$$ is Newton's gravitational constant.
- $$\mu$$ is the string tension, measuring the energy per unit length of the string.

The string's relativistic motion introduces **velocity perturbations** in the wake, given by:

$$
\delta v = 4\pi G\mu v \gamma(v),
$$


where $$v$$ is the string's speed and $$\gamma(v)$$ is the relativistic Lorentz factor. These velocity perturbations result in planar overdensities that appear as a wedge in 3D spacetime.

<details>
<summary><i>Python: Constructing Wake Geometry</i></summary>

import numpy as np
Define cosmic string constants
stringTension = 3E-7 # String tension Gμ
stringSpeed = 0.8 * 299792.458 # String speed in km/s (80% of c)
lorentzFactor = 1 / np.sqrt(1 - (stringSpeed / 299792.458) ** 2) def computeWakeGeometry(wakeLength, wakeDepth, wakeDeficitAngle):
"""
Constructs 3D coordinates for a cosmic string wake.
"""
wakeEndPoints = [
[wakeDepth * np.cos(wakeDeficitAngle / 2), wakeDepth * np.sin(wakeDeficitAngle / 2), 0],
[wakeDepth * np.cos(wakeDeficitAngle / 2), -wakeDepth * np.sin(wakeDeficitAngle / 2), 0],
]
wakeProjection = [
[0, 0, wakeLength],
[wakeEndPoints, wakeEndPoints1
, wakeLength],
[wakeEndPoints1
, wakeEndPoints1
1
, wakeLength],
]
return wakeEndPoints, wakeProjection
Example Wake Geometry
wakeLength = 10.0 # Mpc
wakeDepth = 5.0 # Mpc
wakeDeficitAngle = np.pi / 3
wakeGeometry = computeWakeGeometry(wakeLength, wakeDepth, wakeDeficitAngle)

text

#### Visualization of Wake Geometry:
![Wake Geometry](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/WakeGeometry3D.png?raw=true)

</details>

---

### **2. Embedding Wakes in Redshift Space**

#### Theoretical Background
To analyze the wake in **cosmological observations**, we convert its physical coordinates to **redshift space**. The relationship between distance $$d$$ and redshift $$z$$ is governed by:

$$
d(z) = \int_0^z \frac{c \, dz'}{H(z')},
$$


where $$H(z)$$ is the Hubble parameter in ΛCDM cosmology. This conversion is non-linear and depends on the expansion history of the universe.

Similarly, spatial dilation due to the universe's expansion is accounted for by the **scale factor**:

$$
a(z) = \frac{1}{1+z},
$$


which rescales the wake's transverse dimensions.

<details>
<summary><i>Python: Redshift Conversion</i></summary>

from astropy.cosmology import Planck18
from astropy.cosmology import z_at_value
import astropy.units as u
Define redshift-distance functions
def redshiftToDistance(z):
return Planck18.comoving_distance(z).value def distanceToRedshift(d):
return z_at_value(Planck18.comoving_distance, d * u.Mpc)
Test example
distance = redshiftToDistance(1.0)
redshift = distanceToRedshift(distance)

text

#### Redshift-Distance Relationship:
![Redshift-Scaling](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/RedshiftScaling.png?raw=true)

</details>

---

### **3. Signal Extraction via Matched Filtering**

#### Theoretical Background
A **matched filter** determines the presence of a known signal pattern (the wake gradient) within noisy data. For a 1D signal, the filter is defined as:

$$
S(t) = \sum_k h(t-k) \cdot d(k),
$$


where:
- $$h$$ is the signal template (wake gradient).
- $$d$$ is the noisy data array.

For 2D data, the filter generalizes to:

$$
S(x, y) = \iint h(u, v) d(x-u, y-v) \, du \, dv.
$$


#### Python Implementation:
<details>
<summary><i>Python: 1D Match Filtering</i></summary>

Define a matched filter
def oneDimensionalMatchFilter(dataArray, signalArray):
return np.correlate(dataArray, signalArray, mode="full")
Example inputs and convolution
signalArray = np.linspace(0, 1, 100)
noiseArray = np.random.normal(0, 1, 100)
combinedArray = signalArray + noiseArray matchFilterResult = oneDimensionalMatchFilter(combinedArray, signalArray)

text

#### Match Filter Output:
![Match Filtering](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/1DMatchFilter.png?raw=true)

</details>

---

## Caveats

1. **Thin Wake Limitation**: For realistic string tensions ($$G\mu \sim 10^{-7}$$), the wake is extremely thin and requires high-resolution numerical simulations.
2. **Scaling Instabilities**: Converting physical coordinates to redshift coordinates introduces numerical inaccuracies at extreme $$z$$ values.
3. **Noise Properties**: The current noise models are Gaussian and do not include non-linear perturbations, limiting their realism.

---

## Next Steps

- [ ] Extend the framework to account for **higher-dimensional topological defects** (e.g., textures).  
- [ ] Expand noise simulations to include **non-linear large-scale perturbations**.
- [ ] Optimize signal extraction with **wavelet transforms** for multi-scale analysis.  

---

## References

1. Kibble, T.W.B., "Phase Transitions in the Early Universe," *Journal of Physics A*, 1985.  
2. Brandenberger, R., "Topological Defects and Anisotropies," *Modern Physics Letters A*, 2020.  
3. Morais, J., Blamart, M., "Signal Extraction Methods for Cosmic Strings," 2024.  
