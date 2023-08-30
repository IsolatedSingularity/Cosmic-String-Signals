# Topological Defect Signals
###### Under the supervision of [Professor Robert Brandenberger](https://www.physics.mcgill.ca/~rhb/) at [McGill University](https://www.mcgill.ca/). Funded by the [NSERC USRA](https://www.nserc-crsng.gc.ca/students-etudiants/ug-pc/usra-brpc_eng.asp) award with [FRQNT supplements](https://frq.gouv.qc.ca/en/program/supplements-of-the-nserc-undergraduate-student-research-awards-usra-bpca-2023-2024/).

## Objective

Phase transitions in the very early universe cause spontaneous symmetry breaking events which create topological defects. Linear topological defects associated with the U(1) symmetry group are known as **cosmic strings** and create high energy signals in the form of wakes as they propagate through spacetime. Their presence thus is of interest to study the high energy structure of the universe and the standard model.

The strings occur in a class of renormalizable quantum field theories and are stable under specific conditions on the spacetime manifold's homotopy groups. Their dynamics are given by the Nambu-Gotto action, much like bosonic excitations in string theory. What's more is that gravitational backreaction during the early universe causes primordial ΛCDM noise which hides the string signal. Thus the purpose of this repository is to develop statistics to efficiently extract the cosmic string signal admist the non-linear noise through the framework of 21cm cosmology.

## Code Functionality

*Cosmic String Extraction Statistics.py* builds the cosmic string signal from scratch as a finite density of energy radiating a certain temperature difference admist the 21cm background temperature map. This propagates through spacetime and traces out a wake which has the temperature gradient defined on its convex hull. Then using universe simulations done with 21cmFAST, the string signal is embedded in the primordial noise. Finally, to extract the dynamic signal we make use of statistics such as correlation functions, matched filters, and wavelets. The output are plots of these statistics when the signal of the string is detected in the noise.

The report covering all the theory and code can be found in the main repository as a PDF file.

## Limitations

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

At this point, the problem statement has been defined and the foundational code base has been laid. Potential additions to the code would include higher dimensional topological defect signals, and including an algorithm to invert the redshift function without the limitation of convergence.

![alt text](https://github.com/IsolatedSingularity/Cosmic-String-Wakes/blob/main/Plots/2DConvolution.png?raw=true)
