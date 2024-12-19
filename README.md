# Topological Defect Signals

###### Research project under the supervision of [Professor Robert Brandenberger](https://www.physics.mcgill.ca/~rhb/) at [McGill University](https://www.mcgill.ca/). Collaboration with [Mattéo Blamart](https://inspirehep.net/authors/2077637). Funded by the [NSERC USRA](https://www.nserc-crsng.gc.ca/students-etudiants/ug-pc/usra-brpc_eng.asp) award and [FRQNT supplements](https://frq.gouv.qc.ca).

---

### Objective

Phase transitions in the early universe induced **spontaneous symmetry breaking**, forming **topological defects**. Among these, **cosmic strings**—predicted by $U(1)$ symmetry breaking—leave detectable temperature gradients in spacetime due to the wakes they trace. 

However, their detection is hindered by ΛCDM primordial noise in 21cm cosmology. The aim of this repository is to **construct statistical tools to extract cosmic string signals** hidden amid this background noise. The tools include convex geometry modeling, match filtering, and correlation functions to analyze 21cm redshift data.

---

## Code Functionality

### **1. Cosmic String Wake Geometry in Physical Space**

A cosmic string propagating through spacetime generates a wedge-shaped **wake**. The properties of this wake are determined by:
- Velocity of the string $$v$$.
- String tension $$G\mu$$, where $$\alpha = 8\pi G\mu$$ is the **deficit angle**.

#### Key Equations:
1. **Deficit angle**:  
   $$\alpha = 8\pi G\mu$$  
   $$\delta v = 4\pi G\mu v \gamma(v)$$

2. **Brightness Temperature**:  
   $$  
   \delta T_b(z) = A \cdot \left(1 - \frac{T_{\gamma}(z)}{T_k(z)}\right) \cdot \sqrt{1+z}
   $$

where:
- $$T_{\gamma}(z)$$ is the CMB photon temperature.
- $$T_k(z)$$ is the kinetic temperature of baryons in the wake.

<details>
<summary><i>Python: Wake Geometry Construction</i></summary>

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
Defining wake geometry
def cosmic_string_wake(length, depth, tip, angle):
endpoints = [
tip + depth * np.array([np.cos(angle / 2), np.sin(angle / 2), 0]),
tip + depth * np.array([np.cos(angle / 2), -np.sin(angle / 2), 0]),
]
return np.vstack([tip, *endpoints])
Parameters for visualization
wake_length = 10 # Mpc
wake_depth = 5 # Mpc
wake_tip = np.array([2, 2, 0])
string_deficit_angle = np.pi / 3 wake_wedge = cosmic_string_wake(wake_length, wake_depth, wake_tip, string_deficit_angle)
Convex Hull for visualization
hull = ConvexHull(wake_wedge)
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection='3d')
ax.scatter(wake_wedge[:, 0], wake_wedge[:, 1], wake_wedge[:, 2]) for simplex in hull.simplices:
ax.plot(wake_wedge[simplex, 0], wake_wedge[simplex, 1], wake_wedge[simplex, 2], 'r-') ax.set_title("Cosmic String Wake Geometry")
plt.show()

text

#### Visualization of Wake Geometry:
![Wake Geometry](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/WakeGeometry3D.png?raw=true)

</details>

---

### **2. Embedding Wake in Redshift Space**

To observe signals in cosmology, the wake geometry must be transferred to **redshift space**. Redshift and distance are related through cosmology-dependent scaling relations:
$$
d(z) = \int_0^z \frac{c \, dz'}{H(z')}, \quad a(z) = \frac{1}{1+z}
$$


<details>
<summary><i>Python: Redshift Embedding</i></summary>

from astropy.cosmology import Planck18 as cosmo
Compute redshift-distance relationship
z_vals = np.linspace(0, 4, 500)
distances = cosmo.comoving_distance(z_vals).value plt.figure(figsize=(10, 5))
plt.plot(z_vals, distances)
plt.title("Redshift-Distance Relationship")
plt.xlabel("Redshift (z)")
plt.ylabel("Distance (Mpc)")
plt.grid()
plt.show()

text

#### Redshift-Distance Relationship:
![Redshift to Distance](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/RedshiftScaling.png?raw=true)

</details>

---

### **3. Signal Processing with Match Filtering**

The cosmic string wake generates a narrow **temperature gradient** signal, which we extract using **match filtering**. Match filtering compares a signal template $$h$$ against data $$d$$ for a best match:
$$
S(t) = \sum_k h(t-k) \cdot d(k)
$$


<details>
<summary><i>Python: Match Filtering on Wake Signal</i></summary>

def match_filter(data, template):
return np.correlate(data, template, mode='full')
Example signal and noisy data
wake_signal = np.linspace(0, 1, 100)
noise = np.random.normal(0, 0.1, 100)
data = wake_signal + noise
Match filtering
filter_result = match_filter(data, wake_signal) plt.figure(figsize=(8, 6))
plt.plot(filter_result, label="Match Filter Output")
plt.legend()
plt.title("Signal Extraction with Match Filtering")
plt.show()

text

#### Match Filtering:
![Match Filtering](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/1DMatchFilter.png?raw=true)

</details>

---

### **4. Combining Signal with ΛCDM Noise**

Cosmological ΛCDM noise is simulated using **21cmFAST**. The wake signal is superimposed on a 3D primordial noise map. 

#### Visualization:
![Combined Maps](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/CombinedMapSlice.png?raw=true)

---

## Caveats

1. **Wake Planarity**: For realistic $$G\mu$$ values, the wake becomes nearly planar, which complicates convex hull detection.
2. **Redshift Scaling Instability**: Numerical scaling functions may diverge for very small or very large $$z$$.

---

## Next Steps

- [ ] Extend framework to higher-dimensional defects (e.g., textures).
- [ ] Build GPU-accelerated match filtering for real-time signal processing.
- [ ] Address numerical instability in redshift scaling for extreme $$z$$ values.

For theoretical derivations and complete methods, refer to the full report: [Signal Extraction of Cosmic String Topological Defects](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Signal-extraction-of-cosmic-string-topological-defects-Morais.pdf).
