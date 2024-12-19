# Topological Defect Signals
###### Under the supervision of [Professor Robert Brandenberger](https://www.physics.mcgill.ca/~rhb/) at [McGill University](https://www.mcgill.ca/), and in collaboration with [Mattéo Blamart](https://inspirehep.net/authors/2077637). Funded by the [NSERC USRA](https://www.nserc-crsng.gc.ca/students-etudiants/ug-pc/usra-brpc_eng.asp) award with [FRQNT supplements](https://frq.gouv.qc.ca/en/program/supplements-of-the-nserc-undergraduate-student-research-awards-usra-bpca-2023-2024/).

![alt text](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/2DConvolution.png?raw=true)

## Objective
# Topological Defect Signals
###### Under the supervision of [Professor Robert Brandenberger](https://www.physics.mcgill.ca/~rhb/) at [McGill University](https://www.mcgill.ca/), and in collaboration with [Mattéo Blamart](https://inspirehep.net/authors/2077637). Funded by the [NSERC USRA](https://www.nserc-crsng.gc.ca/students-etudiants/ug-pc/usra-brpc_eng.asp) award with [FRQNT supplements](https://frq.gouv.qc.ca/en/program/supplements-of-the-nserc-undergraduate-student-research-awards-usra-bpca-2023-2024/).

![alt text](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/2DConvolution.png?raw=true)

## Objective

Phase transitions in the very early universe cause spontaneous symmetry breaking events which create topological defects. Linear topological defects associated with the U(1) symmetry group are known as **cosmic strings** and create high energy signals in the form of wakes as they propagate through spacetime. Their presence thus is of interest to study the high energy structure of the universe and the standard model.

The strings occur in a class of renormalizable quantum field theories and are stable under specific conditions on the spacetime manifold's homotopy groups. Their dynamics are given by the Nambu-Gotto action, much like bosonic excitations in string theory. What's more is that gravitational backreaction during the early universe causes primordial ΛCDM noise which hides the string signal. Thus the purpose of this repository is to develop statistics to efficiently extract the cosmic string signal admist the non-linear noise through the framework of 21cm cosmology.

## Code Functionality

*Cosmic String Extraction Statistics.py* builds the cosmic string signal from scratch as a finite density of energy radiating a certain temperature difference admist the 21cm background temperature map. This propagates through spacetime and traces out a wake which has the temperature gradient defined on its convex hull. Then using universe simulations done with 21cmFAST, the string signal is embedded in the primordial noise. Finally, to extract the dynamic signal we make use of statistics such as correlation functions, matched filters, and wavelets. The output are plots of these statistics when the signal of the string is detected in the noise.

The report covering all the theory and code can be found in the main repository as a PDF file.

## Caveats

The method in which points are detected within the wake is done using complex convex hulls. This algorithm
becomes problematic when the blown up deficit angle is replaced by its actual value of $\alpha = 8 \pi G \mu$ which
is very small and thus the wake becomes a plane. The algorithm is based on connecting simplices along
different vertices and does not work when the topology of the object is in 1D. Next, when converting from
physical to comoving coordinates, one uses an inverse scaling factor of the form $a^{−1}(z) = (1 − z)/z$, which
can also be substituted for $a^{−1}(t_0) \sim 10^3$ for current observations. This scaling becomes an issue when
wanting to scale physical axes to redshift axes using the numerical function from the astropy package, which
doesn’t converge for small $\mathcal{O}(1)$ or large $\mathcal{O}(1000)$ values of redshift. Thus, we are left with to work in a
snapshot of physical coordinates to substitute for a continuous comoving coordinate system.

## Next Steps

At this stage, we've clarified the problem statement and established the fundamental code framework. Potential additions to the code would include higher dimensional topological defect signals, and including an algorithm to invert the redshift function without the limitation of convergence.
Phase transitions in the very early universe cause spontaneous symmetry breaking events which create topological defects. Linear topological defects associated with the U(1) symmetry group are known as **cosmic strings** and create high energy signals in the form of wakes as they propagate through spacetime. Their presence thus is of interest to study the# Topological Defect Signals
###### Under the supervision of [Professor Robert Brandenberger](https://www.physics.mcgill.ca/~rhb/) at [McGill University](https://www.mcgill.ca/), and in collaboration with [Mattéo Blamart](https://inspirehep.net/authors/2077637). Funded by the [NSERC USRA](https://www.nserc-crsng.gc.ca/students-etudiants/ug-pc/usra-brpc_eng.asp) award with [FRQNT supplements](https://frq.gouv.qc.ca/en/program/supplements-of-the-nserc-undergraduate-student-research-awards-usra-bpca-2023-2024/).

![alt text](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/2DConvolution.png?raw=true)

## Objective

Phase transitions in the very early universe induce spontaneous symmetry breaking, producing topological defects. Linear topological defects associated with the U(1) symmetry group are known as **cosmic strings**, which form during inflationary phase transitions and persist into the late universe. These defects leave long-lived signals in the form of **wakes** as they propagate through spacetime.

Their observational signature manifests as anisotropies in the Cosmic Microwave Background (CMB), particularly detectable through **21 cm cosmology**. However, these signals are hidden amidst complex non-linear primordial noise resulting from $\Lambda$CDM cosmological perturbations. Extracting these signals efficiently is crucial for probing the high-energy structure of the early universe and constraining grand unified theories.

Mathematically, cosmic strings arise from the non-triviality of homotopy groups of the vacuum manifold, typically $\pi_1(M) \neq \identity$, where $M$ is the vacuum manifold of a spontaneously broken symmetry group $G/H$. Their stability conditions can be derived from exact homotopy sequences, and their gravitational lensing leads to planar overdensities known as wakes. These wakes generate a brightness temperature contrast $\delta T_b(z)$ which can be observed in the 21 cm line, allowing us to probe energy scales up to $G \mu \sim 10^{-7}$.

## Code Functionality

Given a 3D spatial snapshot of a Hubble volume at recombination, we construct a cosmic string wake wedge based on a tension parameter $G \mu$, randomize its orientation, and map physical coordinates to redshift space. This wedge is immersed into a $\Lambda$CDM noise field generated by 21cmFAST simulations. By assigning a brightness temperature $\delta T_b$ to points inside the wedge, and superimposing them with primordial noise, we obtain a 3D temperature field.

To isolate the cosmic string signal from noise, we implement **match filtering** techniques. Match filters, defined as convolutional statistics, maximize the signal-to-noise ratio. We unfold 2D slices of the 3D cube along different orientations (horizontal/vertical), apply match filters, and identify the domain in which the cosmic string signal resides.

### Mathematical Formulation

A cosmic string wake is characterized by a deficit angle $\alpha = 8 \pi G \mu$ in the transverse plane. Its geometry in physical space is defined by a wedge with a length scale given by the Hubble length $l_H$ and thickness related to the string formation time $t_i$:

$$
\text{Wake dimensions: } \quad \text{length} \sim c_1 t_i, \quad \text{depth} \sim \gamma_s v_s t_i, \quad \text{width} \sim \alpha \times (\text{some scale})
$$

Here $v_s$ is the string velocity and $\gamma_s$ the corresponding Lorentz factor. Converting distances $d$ to redshifts $z$ uses the Planck18 cosmology:

$$
d(z) = \int_0^z \frac{c \, dz'}{H(z')}, \quad \text{and} \quad H(z) = H_0 \sqrt{\Omega_m(1+z)^3 + \Omega_\Lambda},
$$

enabling us to scale the vertical axis from Mpc to redshift space. The brightness temperature perturbation $\delta T_b(z)$ inside the wake can be approximated by:

$$
\delta T_b(z) = [0.07 K] \frac{x_c}{1+x_c} \left( 1 - \frac{T_\gamma(z)}{T_{K/g}(z)} \right)\sqrt{1+z} \, ,
$$

where $x_c$ is the collision coefficient, $T_\gamma(z)$ is the CMB photon temperature, and $T_{K/g}(z)$ the kinetic/gas temperature inside the wake.

### Code Snippets

Below is a snippet from the main code *Cosmic String Extraction Statistics.py*, showcasing how we construct the wake in physical space and then rotate it using $SO(3)$ transformations.

<details>
  <summary><i>Click to expand Python code snippet</i></summary>

```python
# Defining group representations for SO(3) rotations
xAxisRotation = lambda theta : np.array([
        [1,0,0],
        [0,np.cos(theta),-np.sin(theta)],
        [0,np.sin(theta),np.cos(theta)]
    ])

yAxisRotation = lambda theta : np.array([
        [np.cos(theta),0,np.sin(theta)],
        [0,1,0],
        [-np.sin(theta),0,np.cos(theta)]
    ])

zAxisRotation = lambda theta : np.array([
        [np.cos(theta),-np.sin(theta),0],
        [np.sin(theta),np.cos(theta),-np.sin(theta)],
        [0,0,1]
    ])

# Shifting the wake to the origin, applying random rotations, and shifting it back
shiftedWedge = wakeWedge - wakeMiddlePoint
modifiedWedgeI = np.dot(shiftedWedge, xAxisRotation(np.random.uniform(0,2*np.pi)).T)
modifiedWedgeII = np.dot(modifiedWedgeI, yAxisRotation(np.random.uniform(0,2*np.pi)).T)
modifiedWedgeIII = np.dot(modifiedWedgeII, zAxisRotation(np.random.uniform(0,2*np.pi)).T)
rotatedWedge = modifiedWedgeIII + wakeMiddlePoint

</details>

Following the geometry construction, we apply a convex hull method to determine which points lie inside the wedge and assign each of these points a temperature value from $\delta T_b(z)$.
<details> <summary><i>Click to expand Python code snippet</i></summary>

from scipy.spatial import ConvexHull

# Determining points inside the convex hull
hull = ConvexHull(wakeWedge)
# For a given point set testPoints, we check if they lie inside the hull
def pointChecker(hullPoints, testPoints):
    hull = ConvexHull(hullPoints)
    newPoints = np.vstack((hullPoints, testPoints))
    newHull = ConvexHull(newPoints)
    return list(hull.vertices) == list(newHull.vertices)

# Assigning temperature gradient to points inside the wedge
def gradientManifold(wedge, resolution):
    # Generate mesh grid, check points, assign δT_b if inside wedge
    # Return arrays of wedge points and temperatures
    ...
    return wedgePoints, wedgeTemperatures

</details>
Caveats

    Convex Hull Limitation: For actual cosmic string parameters, $\alpha$ is extremely small, making the wake nearly planar. The convex hull algorithm, based on connecting simplices, struggles when the object’s topology degenerates into lower dimensions.

    Redshift Conversion Issues: When converting physical coordinates to redshift space using astropy’s cosmology functions, convergence issues arise for extremely small $\mathcal{O}(1)$ or large $\mathcal{O}(1000)$ redshifts. We thus work with snapshot approximations rather than continuous mappings.

Next Steps

Implemented basic framework for constructing cosmic string wakes and embedding them in noise.
Developed match filtering techniques for 1D and 2D slices to isolate the string signal.
Extend method to consider more realistic string tensions and integrate machine learning techniques for automated signal detection.
Improve redshift space conversion algorithms for better accuracy at high/low $z$ regimes.

    Consider other topological defects (monopoles, domain walls, textures) by examining $\pi_n(M)$ for $n>1$ and including their corresponding field configurations.

    [!TIP] Refer to the accompanying PDF report for a more in-depth theoretical explanation, including homotopy groups, stability conditions, and the topological quantum field theory framework behind topological defects.

    [!NOTE] Once algorithms are validated, they can be optimized for parallel computing architectures and possibly integrated into large-scale cosmological simulations, improving the search for cosmic strings in future 21 cm surveys.

 high energy structure of the universe and the standard model.

The strings occur in a class of renormalizable quantum field theories and are stable under specific conditions on the spacetime manifold's homotopy groups. Their dynamics are given by the Nambu-Gotto action, much like bosonic excitations in string theory. What's more is that gravitational backreaction during the early universe causes primordial ΛCDM noise which hides the string signal. Thus the purpose of this repository is to develop statistics to efficiently extract the cosmic string signal admist the non-linear noise through the framework of 21cm cosmology.

## Code Functionality

*Cosmic String Extraction Statistics.py* builds the cosmic string signal from scratch as a finite density of energy radiating a certain temperature difference admist the 21cm background temperature map. This propagates through spacetime and traces out a wake which has the temperature gradient defined on its convex hull. Then using universe simulations done with 21cmFAST, the string signal is embedded in the primordial noise. Finally, to extract the dynamic signal we make use of statistics such as correlation functions, matched filters, and wavelets. The output are plots of these statistics when the signal of the string is detected in the noise.

The report covering all the theory and code can be found in the main repository as a PDF file.

## Caveats

The method in which points are detected within the wake is done using complex convex hulls. This algorithm
becomes problematic when the blown up deficit angle is replaced by its actual value of $\alpha = 8 \pi G \mu$ which
is very small and thus the wake becomes a plane. The algorithm is based on connecting simplices along
different vertices and does not work when the topology of the object is in 1D. Next, when converting from
physical to comoving coordinates, one uses an inverse scaling factor of the form $a^{−1}(z) = (1 − z)/z$, which
can also be substituted for $a^{−1}(t_0) \sim 10^3$ for current observations. This scaling becomes an issue when
wanting to scale physical axes to redshift axes using the numerical function from the astropy package, which
doesn’t converge for small $\mathcal{O}(1)$ or large $\mathcal{O}(1000)$ values of redshift. Thus, we are left with to work in a
snapshot of physical coordinates to substitute for a continuous comoving coordinate system.

## Next Steps

At this stage, we've clarified the problem statement and established the fundamental code framework. Potential additions to the code would include higher dimensional topological defect signals, and including an algorithm to invert the redshift function without the limitation of convergence.
