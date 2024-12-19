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
  <summary><i>Cosmic String Signal Extraction Algorithm</i></summary>

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

# Function that quickly plots 3D figures in physical space
def quickPlot(
    scale,points,labels,colors,colormaps,convexhull,testPoints=[]
    ):
    
    # Defining some plotting parameters 
    fig = plt.figure(figsize=(12, 10), dpi=80)
    ax = fig.add_subplot(111, projection="3d")
    
    # Finding the approximate center point of wake to scale the plots properly
    approximateCenter = []
    for i in points:
        if len(i) == 6:
            for j in range(0,3):
                approximateCenter.append(
                    0.5*(np.max(i.T[j])-np.min(i.T[j])) + np.min(i.T[j])
                )

    # Defining axes labels and limits
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('y [Mpc]')
    ax.set_zlabel('z [Mpc]')
    ax.set_xlim(approximateCenter[0]-scale/2.7,approximateCenter[0]+scale/2.7)
    ax.set_ylim(approximateCenter[1]-scale/2.7,approximateCenter[1]+scale/2.7)
    ax.set_zlim(approximateCenter[2]-scale/2.7,approximateCenter[2]+scale/2.7)
    
    # Testing wether or not there are convex hulls involved
    if convexhull == False:
    
        # Differentiating if points have an optional temperature color map
        if len(colormaps) != 0:
            
            physicalPlot = ax.scatter(
                points[-1][:,0],points[-1][:,1],points[-1][:,2], c = colormaps,
                cmap = plt.cm.viridis, s = 90
            )
        
            fig.colorbar(physicalPlot, ax=ax, label = 'δT(z)', pad=0.1)
        
        if len(colormaps) == 0:
            
            # Iterating through sets of points being plotted
            for pointType in range(len(points)):
                
                # Differentiating between individual points and sets of points
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
        
        # Defining hull & plotting hull vertices
        hull = ConvexHull(points[0])
        ax.scatter(
            points[0].T[0], points[0].T[1], points[0].T[2], color='red', s = 90
            )
        
        # Plotting hull simplices that are connected to vertices
        for s in hull.simplices:
            s = np.append(s, s[0]) #simplice
            ax.plot(points[0][s, 0], points[0][s, 1], points[0][s, 2], color='maroon')
        
        # Plotting any potential test points along side the hull
        if len(testPoints) != 0:
            ax.scatter(
            testPoints[:,0],testPoints[:,1],testPoints[:,2], color='black', s = 90
            )
    
    # Design parameters
    ax.legend()
    plt.tight_layout()
    plt.autoscale()
    plt.show()

    return

# Function that quickly plots 3D figures in redshift-physical space
def quickPlotRedshift(
    scale,points,labels,colors,colormaps,convexhull,testPoints=[]
    ):
    
    # Defining some plotting parameters 
    fig = plt.figure(figsize=(12, 10), dpi=80)
    ax = fig.add_subplot(111, projection="3d")
    redshiftScale = 6.77E-5
    
    # Finding the approximate center point of wake to scale the plots properly
    approximateCenter = []
    for i in points:
        if len(i) == 6:
            for j in range(0,3):
                approximateCenter.append(
                    0.5*(np.max(i.T[j])-np.min(i.T[j])) + np.min(i.T[j])
                )

    # Defining axes labels and limits
    ax.set_xlabel('x [Mpc]')
    ax.set_ylabel('y [Mpc]')
    ax.set_zlabel('z')
    ax.set_xlim(approximateCenter[0]-scale/2,approximateCenter[0]+scale/2)
    ax.set_ylim(approximateCenter[1]-scale/2,approximateCenter[1]+scale/2)
    ax.set_zlim(approximateCenter[2]-redshiftScale/2,approximateCenter[2]+redshiftScale/2)
    
    # Testing wether or not there are convex hulls involved
    if convexhull == False:
    
        # Differentiating if points have an optional temperature color map
        if len(colormaps) != 0:
            
            physicalPlot = ax.scatter(
                points[-1][:,0],points[-1][:,1],points[-1][:,2], c = colormaps,
                cmap = plt.cm.viridis, s = 90
            )
        
            fig.colorbar(physicalPlot, ax=ax, label = 'δT(z)', pad=0.1)
        
        if len(colormaps) == 0:
            
            # Iterating through sets of points being plotted
            for pointType in range(len(points)):
                
                # Differentiating between individual points and sets of points
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
        
        # Defining hull & plotting hull vertices
        hull = ConvexHull(points[0])
        ax.scatter(
            points[0].T[0], points[0].T[1], points[0].T[2], color='red', s = 90
            )
        
        # Plotting hull simplices that are connected to vertices
        for s in hull.simplices:
            s = np.append(s, s[0]) #simplice
            ax.plot(points[0][s, 0], points[0][s, 1], points[0][s, 2], color='maroon')
        
        # Plotting any potential test points along side the hull
        if len(testPoints) != 0:
            ax.scatter(
            testPoints[:,0],testPoints[:,1],testPoints[:,2], color='black', s = 90
            )
    
    # Design parameters
    ax.legend()
    plt.tight_layout()
    plt.autoscale()
    plt.show()

    return


# %% Defining the Hubble volume and cosmic string wake wedge 유

# Initializing required parameters for the universe model
hubbleParameter = 70.41*1000 #[m/s/Mpc] H(z=0)=H0
densityRatio = 0.3 #non-relativistic matter density for flat ΛCDM 
speedOfLight = 299792458 #[m/s]
hubbleScale = speedOfLight/hubbleParameter #[Mpc] Hubble length
hubbleLattice = hubbleScale * np.array([
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
     ])

# Defining required parameters for the string model
Gμ = 3E-7 #string tension
deficitAngle = np.pi/3 #actually is 8*np.pi*Gμ, using π/3 to visualize the angle
c1 = rd.uniform(1-0.01,1+0.01) #O(1) constant
meterToMpc = 3.24078E-23 #1m = 3.24078E-23 Mpc; to convert m -> Mpc or rather m/Mpc -> 1 
redshiftTime = lambda z : 2/(3 * hubbleParameter * (densityRatio**(1/2)) * (1+z)**(3/2)) * meterToMpc**(-1) #t(z) in s
formationTime = redshiftTime(1100) #time during which string formation is most prominent: recombination z~1100

# Defining the wake geometry in 3D physical space
stringSpeed = speedOfLight * rd.uniform(0.5,0.9) #[m/s] conservative estimate (50%-90% of c)
gammaFactor = ( 1/( 1 - (stringSpeed/speedOfLight)**2 ) )**(1/2)
wakeLength = c1 * formationTime * speedOfLight * meterToMpc #[Mpc]
wakeDepth = formationTime * gammaFactor * stringSpeed * meterToMpc #[Mpc] radial length
wakeTipPoint = np.array([
    rd.uniform(0,0.9*hubbleScale),rd.uniform(0,0.9*hubbleScale),rd.uniform(0,0.9*hubbleScale)
])
wakeEndPoints = np.array([
    [wakeTipPoint[0]+wakeDepth*np.cos(deficitAngle/2),wakeTipPoint[1]+wakeDepth*np.sin(deficitAngle/2),wakeTipPoint[2]],
    [wakeTipPoint[0]+wakeDepth*np.cos(deficitAngle/2),wakeTipPoint[1]-wakeDepth*np.sin(deficitAngle/2),wakeTipPoint[2]]
    ])
projectedWakePoints = np.array([
    wakeTipPoint+[0,0,wakeLength],wakeEndPoints[0]+[0,0,wakeLength],wakeEndPoints[1]+[0,0,wakeLength]
    ])
wakeMiddlePoint = np.array([
    wakeTipPoint[0]+(1/2)*wakeDepth*np.cos(deficitAngle/2),wakeTipPoint[1],wakeTipPoint[2]+wakeLength/2
    ])
wakeWedge = np.array([
    wakeTipPoint, wakeEndPoints[0], wakeEndPoints[1], projectedWakePoints[0], projectedWakePoints[1], projectedWakePoints[2]
])

scaleConverter = hubbleScale/wakeLength * 1/2 #used to scale plot sizes

# Plotting Hubble lattice and wake wedge
quickPlot(
    hubbleScale,[hubbleLattice, wakeWedge],['Hubble lattice','Wake wedge'],
    ['blue','red'],colormaps=[],convexhull=False
    )

# Plotting again with smaller axis
quickPlot(
    hubbleScale/scaleConverter,[wakeWedge, wakeMiddlePoint],['Wake wedge','Wake center'],
    ['red','black'],colormaps=[],convexhull=False
    )


# %% Generating rotations via group representations and translations on wake wedge state 유

# Real representation of SO(3) group
xAxisRotation = lambda theta : np.array([
        [1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]
    ])

yAxisRotation = lambda theta : np.array([
        [np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]
    ])

zAxisRotation = lambda theta : np.array([
        [np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),-np.sin(theta)],[0,0,1]
    ])

# Shifting the wake to the origin to then perform rotations
shiftedWedge = np.add(wakeWedge,-wakeMiddlePoint)

# Applying rotations on wake (random angles around all axes)
modifiedWedgeI = np.dot(shiftedWedge, xAxisRotation(rd.uniform(0,2*np.pi)).T)
modifiedWedgeII = np.dot(modifiedWedgeI, yAxisRotation(rd.uniform(0,2*np.pi)).T)
modifiedWedgeIII = np.dot(modifiedWedgeII, zAxisRotation(rd.uniform(0,2*np.pi)).T)

# Shifting wake back to center created
rotatedWedge = np.add(modifiedWedgeIII,wakeMiddlePoint)

# Plotting unrotated and rotated wake in the Hubble lattice
quickPlot(
    hubbleScale/scaleConverter,[rotatedWedge, wakeWedge, wakeMiddlePoint],['Rotated wake','Wake wedge','Wake center'],
    ['darkorchid','red','black'],colormaps=[],convexhull=False
    )


#%% Computations to convert z physical coordinates in Mpc to redshift axis 유

# Distance in Mpc as a function of redshift & its inverse
redshiftToDistance = lambda z : Planck18.comoving_distance(z).value #[Mpc] d(z)
distanceToRedshift = lambda d : z_at_value(Planck18.comoving_distance,(d*u.Mpc)).value # z(d)

# Plotting the relation functions
x = np.linspace(5,100,100) #redshift points
y = np.linspace(5,12000,100) #physical points in Mpc

fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

# Distance as a function of redshift plot
ax = axs[0]
ax.plot(x, redshiftToDistance(x), color='blue', label='d(z)')
plt.setp(ax, xlabel='z' , ylabel='Distance [Mpc]')
ax.legend()
ax.grid()

# Redshift as a function of distance plot
ax = axs[1]
ax.plot(y, distanceToRedshift(y), color='red', label='z(d)')
plt.setp(ax, ylabel='z' , xlabel='Distance [Mpc]')
ax.grid()
ax.legend()

# Demonstrating inverse function equivalence
print('d(z=100): ', int(redshiftToDistance(100)), 'Mpc')
print('z(d~12831 Mpc): ', int(distanceToRedshift(12831.379081434678)))

# Defining z(d) scaling function to apply to arbitrary wake wedges
def wedgeToRedshift(wedge):
    
    wedgeCopy = np.copy(wedge) #to not replace original input wedge
    
    # Checking if using one point or a set of points
    if len(wedge) == 3: #point
        wedgeCopy[-1] = distanceToRedshift(wedge[-1])
        
    else:
        
        for point in range(len(wedge)):
            wedgeCopy[point][-1] = distanceToRedshift(wedge[point][-1])
            
    return wedgeCopy

scaledWedge = wedgeToRedshift(wakeWedge)
scaledRotatedWedge = wedgeToRedshift(rotatedWedge)
scaledWedgeCenter = wedgeToRedshift(wakeMiddlePoint)

# Plotting scaled (d [Mpc] -> z) wedges in 2D physical x 1D redshift space aka (2,1) space
quickPlotRedshift(
    hubbleScale/scaleConverter,[scaledRotatedWedge, scaledWedge, scaledWedgeCenter],['Rotated wake','Wake wedge','Wake center'],
    ['darkorchid','red','black'],colormaps=[],convexhull=False
    )


#%% Computations to convert comoving coordinates to expanding physical coordinates 유

scaleFactor = lambda z : 1/(1+z) #a(t) = 1/(1+z)

# Function that scales xy plane via the scale factor
def comovingToPhysical(wedge):
    
    wedgeCopy = np.copy(wedge)
    
    # Checking if using one point or a set of points
    if len(wedge) == 3: #point
        wedgeCopy[0] = scaleFactor(wedge[-1])*wedgeCopy[0] #x_phys = a(t)*x_com
        wedgeCopy[1] = scaleFactor(wedge[-1])*wedgeCopy[1]
        
    else:
        
        for point in range(len(wedge)):
            wedgeCopy[point][0] = scaleFactor(wedge[point][-1])*wedgeCopy[point][0]
            wedgeCopy[point][1] = scaleFactor(wedge[point][-1])*wedgeCopy[point][1]
            
    return wedgeCopy

physicalScaledWedge = comovingToPhysical(scaledWedge)
physicalScaledRotatedWedge = comovingToPhysical(scaledRotatedWedge)
physicalScaledWedgeCenter = comovingToPhysical(scaledWedgeCenter)

# Plotting scaled (d [Mpc] -> z) wedges in 2D physical x 1D redshift space aka (2,1) space
quickPlotRedshift(
    hubbleScale/scaleConverter,[physicalScaledRotatedWedge, physicalScaledWedge, physicalScaledWedgeCenter],
    ['Physical rotated wake','Physical Wake wedge','Physical wake center'],
    ['darkorchid','red','black'],colormaps=[],convexhull=False
    )


#%% Computing the boundary of the wedge via convex hulls & calculating if test points are within it 유

# Plotting convex hulls for different variants of the wedge in space
quickPlot(
    hubbleScale/scaleConverter,[wakeWedge,wakeWedge],[],[],colormaps=False,convexhull=True
    )
quickPlotRedshift(
    hubbleScale/scaleConverter,[physicalScaledWedge,physicalScaledWedge],[],[],colormaps=False,convexhull=True
    )
quickPlotRedshift(
    hubbleScale/scaleConverter,[physicalScaledRotatedWedge,physicalScaledRotatedWedge],[],[],colormaps=False,convexhull=True
    )

# Function that checks if a set of points lie within the hull
def pointChecker(scale,hullPoints, testPoints):
    
    # Defining hull and test points
    hull = ConvexHull(hullPoints)
    newPoints = np.append(hullPoints, testPoints, axis=0)
    newHull = ConvexHull(newPoints)
    
    # Checking if points are within the hull
    if list(hull.vertices) == list(newHull.vertices):
        criterion = True
    else:
        criterion = False
        
    return criterion
    
# Running preliminary examples
pointChecker(hubbleScale/scaleConverter,wakeWedge,np.array([[0.1,0.15,0.05]]))
print(pointChecker(hubbleScale/scaleConverter,wakeWedge,np.array([[0.1,0.15,0.05]])))

pointChecker(hubbleScale/scaleConverter,physicalScaledRotatedWedge,np.array([[0.2,0.1,1E-5]]))
print(pointChecker(hubbleScale/scaleConverter,physicalScaledRotatedWedge,np.array([[0.2,0.1,1E-5]])))


#%% Computing brightness temperature function δT(z) for points within wake 유

# Defining functions for brightness temperature function
photonCMBTemp = lambda z : 2.72548*(1+z) #[K]
kineticTemperature = lambda z : 20*(Gμ*10**6)**2*(gammaFactor*stringSpeed)**2*(1100+1)/(z+1) #[K]
gasTemperature = lambda z : 0.02*(z+1)**2

def brightnessTemperature(z):
    deexcitationCrossSection = 0.16
    if (kineticTemperature(z)>3*gasTemperature(z)):
         u = kineticTemperature(z)
    else:
        u = 3*gasTemperature(z)
    if (kineticTemperature(z)>3*gasTemperature(z)):
        c=4
    else:
        c = 1+kineticTemperature(z)/gasTemperature(z)
    a = 0.017*deexcitationCrossSection/(1+deexcitationCrossSection)*(1-photonCMBTemp(z)/u)*np.sqrt(1+z)*c
    return a

brightnessDistanceTemperature = lambda d : brightnessTemperature(distanceToRedshift(d)) #[K] δT(d)

# Old brightness temperature function definitions
# kineticTemperature = lambda z : (22020*Gμ*(1E12)*(stringSpeed*gammaFactor)**2*speedOfLight**(-2))/(z+1) #[K]
# deexcitationCrossSection = lambda z : 2.76E-9 * (1 - 1/((1+0.0007112 * (kineticTemperature(z))**(2.28) )**(0.014)))*(1/(100**3))
# photonCMBTemp = lambda z : 2.72548*(1+z) #[K] from wiki
# densityNumber = lambda z : 0.76* (1+z)**3 #[1/m^3]
# collisionCoefficient = lambda z : (densityNumber(z)*deexcitationCrossSection(z)*0.068)/(2.80E-15 * photonCMBTemp(z)) #[0]
# brightnessTemperature = lambda z : 0.07 * (collisionCoefficient(z)/(1+collisionCoefficient(z))) * (1-photonCMBTemp(z)/kineticTemperature(z)) * np.sqrt(1+z) #[K] δT(z)
# brightnessDistanceTemperature = lambda d : brightnessTemperature(distanceToRedshift(d)) #[K] δT(d)

# Defining function to find points within wedge and apply temperature gradient to those points
def gradientManifold(scale,resolution,wedge,redshiftMpc=[]):
    
    # Finding the approximate center point of wake to scale the plots properly
    approximateCenter = []
    redshiftScale = 6.77E-5
    
    if len(wedge) == 6:
        for j in range(0,3):
            approximateCenter.append(
                np.mean(wedge.T[j])
            )
            
    # Determining if z axis is scaled to redshift, and if so what physicalScale should be used for the xy plane
    if len(redshiftMpc) == 0: #Mpc scale, not redshift

        # Finding grid of points within the Hubble lattice
        meshGrid = np.mgrid[
            approximateCenter[0]-scale/2:approximateCenter[0]+scale/2:resolution, 
            approximateCenter[1]-scale/2:approximateCenter[1]+scale/2:resolution, 
            approximateCenter[2]-scale/2:approximateCenter[2]+scale/2:resolution
            ]
        X, Y, Z = meshGrid
        positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
    
    else:
        
        # Finding grid of points within the Hubble lattice
        meshGrid = np.mgrid[
            approximateCenter[0]-scale/2:approximateCenter[0]+scale/2:resolution, 
            approximateCenter[1]-scale/2:approximateCenter[1]+scale/2:resolution, 
            approximateCenter[2]-redshiftScale/2:approximateCenter[2]+redshiftScale/2:1.67E-6
            ]
        X, Y, Z = meshGrid
        positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
        
    # Getting array of all points within the space
    emptyPositions = []
    wedgePoints = [] #Physical points within wedge only
    wedgeTemperatures = [] #Temperature of points within wedge only

    # Converting mesh grid to physical points in 3D space
    for i in range(len(positions.T)):
        emptyPositions.append([positions[0][i],positions[1][i],positions[2][i]])
    physicalPositions = np.array(emptyPositions)

    # Iterating through the Hubble lattice to find points within wake and assign a temperature
    temperatureValues = np.zeros((len(positions.T),1))   
        
    # Assigning a temperature gradient to points within wake wedge    
    wakeCopy = np.copy(wedge)
    
    # Determining if z axis is scaled to redshift   
    if len(redshiftMpc) == 0:
        for i in range(len(physicalPositions)):
            if pointChecker(hubbleScale/scaleConverter,wakeCopy,np.array([physicalPositions[i]])) == True:
                temperatureValues[i] = brightnessDistanceTemperature(physicalPositions[i][2])
                wedgePoints.append(physicalPositions[i])
                wedgeTemperatures.append(brightnessDistanceTemperature(physicalPositions[i][2]))    
            
    else:
        for i in range(len(physicalPositions)):
            if pointChecker(hubbleScale/scaleConverter,wakeCopy,np.array([physicalPositions[i]])) == True:
                temperatureValues[i] = brightnessTemperature(physicalPositions[i][2])
                wedgePoints.append(physicalPositions[i])
                wedgeTemperatures.append(brightnessTemperature(physicalPositions[i][2]))   
        
    physicalWedgePoints = np.array(wedgePoints)
    physicalWedgeTemperatures = np.array(wedgeTemperatures)
    
    return (physicalWedgePoints,physicalWedgeTemperatures,physicalPositions,temperatureValues)

wedgeGradient = gradientManifold(hubbleScale/scaleConverter,0.01,wakeWedge)
rotatedWedgeGradient = gradientManifold(hubbleScale/scaleConverter,0.01,physicalScaledRotatedWedge,redshiftMpc=[9])

# Plotting physical points within wake excluding ambient points
quickPlot(
    hubbleScale/scaleConverter,[wakeWedge,wedgeGradient[0]],['Internal temperature points'],
    ['red'], colormaps=wedgeGradient[1],convexhull=False
    )

quickPlotRedshift(
    hubbleScale/scaleConverter,[physicalScaledRotatedWedge,rotatedWedgeGradient[0]],['Internal temperature points'],
    ['red'], colormaps=rotatedWedgeGradient[1],convexhull=False
    )


#%% Unfolding data in different orientations 유

# Defining & plotting example grids for match filtering
randomValues = np.random.random((60, 60)) #noise map

gradientValues = np.zeros(60*60).reshape(60,60) #wake signal map
for j in range(len(gradientValues)):
    for i in range(len(gradientValues)):
        gradientValues[j][i] = j/60
        
lineValues = np.zeros(60*60).reshape(60,60) #test line signal map
for j in range(len(lineValues)):
    if (j == 29) or (j == 30):
        lineValues.T[j] = 1
    else:
        lineValues.T[j] = 0

# Plotting results of unfolding (without noise)
fig, axs = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.imshow(gradientValues)
ax.set_title('Gradient Map')

ax = axs[1]
ax.imshow(randomValues)
ax.set_title('Noise Map')

ax = axs[2]
ax.imshow(lineValues)
ax.set_title('Line Map')

fig.colorbar(ax.imshow(lineValues), ax=ax, fraction=0.05, pad=0.04, label = 'δT')


# Defining horizontal/vertical unfolding function for 1D match filter
def unfolder(grid,type):
    amplitudeValues = []
    virtualGrid = np.copy(grid)
    
    # Differentiating between unfolding orientation
    if type == 'horizontal':
        for row in virtualGrid:
            for value in row:
                amplitudeValues.append(value)
        return amplitudeValues
    
    if type == 'vertical':
        for column in virtualGrid.T:
            for value in column:
                amplitudeValues.append(value)
        return amplitudeValues

# Unfolding maps in both orientations
horizontalUnfoldedRandom = unfolder(randomValues,'horizontal')
verticalUnfoldedRandom = unfolder(randomValues,'vertical')
horizontalUnfoldedGradient = unfolder(gradientValues, 'horizontal')
verticalUnfoldedGradient = unfolder(gradientValues, 'vertical')
horizontalUnfoldedLine = unfolder(lineValues,'horizontal')
verticalUnfoldedLine = unfolder(lineValues,'vertical')

unfoldSize = len(horizontalUnfoldedGradient)
unfoldedIndices = np.linspace(0,unfoldSize,unfoldSize)
indicesSize = len(unfoldedIndices)/2

# Plotting results of unfolding (without noise)
fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    unfoldedIndices, horizontalUnfoldedGradient, color='darkviolet', label='Gradient'
)
ax.plot(
    unfoldedIndices, horizontalUnfoldedLine, color='mediumspringgreen', label='Line'
)
plt.setp(ax, xlabel='Index' , ylabel='Brightness Temperature [K]')
ax.set_xlim(indicesSize-200,indicesSize+200)
ax.set_title('Horizontal Unfolding')
ax.legend()
ax.grid()

ax = axs[1]
ax.plot(
    unfoldedIndices, verticalUnfoldedGradient, color='darkviolet', label='Gradient'
)
ax.plot(
    unfoldedIndices, verticalUnfoldedLine, color='mediumspringgreen', label='Line'
)
plt.setp(ax, xlabel='Index')
ax.set_xlim(indicesSize-200,indicesSize+200)
ax.set_title('Vertical Unfolding')
ax.legend()
ax.grid()

# Plotting results of unfolding (with noise)
fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    unfoldedIndices, horizontalUnfoldedRandom, color='royalblue', label='Noise'
)
plt.setp(ax, xlabel='Index' , ylabel='Brightness Temperature [K]')
ax.set_xlim(indicesSize-200,indicesSize+200)
ax.set_title('Horizontal Unfolding')
ax.legend()
ax.grid()

ax = axs[1]
ax.plot(
    unfoldedIndices, verticalUnfoldedRandom, color='royalblue', label='Noise'
)
plt.setp(ax, xlabel='Index')
ax.set_xlim(indicesSize-200,indicesSize+200)
ax.set_title('Vertical Unfolding')
ax.legend()
ax.grid()

#%% Match filtering: 1D Examples 유

# Defining (horizontal/vertical)(template)(data) & plotting match filters
horizontalGradientGradient = np.correlate(
    horizontalUnfoldedGradient,horizontalUnfoldedGradient,mode='full'
)
horizontalGradientRandom = np.correlate(
    horizontalUnfoldedGradient,horizontalUnfoldedRandom,mode='full'
)
horizontalGradientLine = np.correlate(
    horizontalUnfoldedGradient,horizontalUnfoldedLine,mode='full'
)
verticalGradientGradient = np.correlate(
    verticalUnfoldedGradient,verticalUnfoldedGradient,mode='full'
)
verticalGradientRandom = np.correlate(
    verticalUnfoldedGradient,verticalUnfoldedRandom,mode='full'
)
verticalGradientLine = np.correlate(
    verticalUnfoldedGradient,verticalUnfoldedLine,mode='full'
)

matchFilterSize = len(horizontalGradientGradient)
filterLinSpace = np.linspace(0,matchFilterSize,matchFilterSize)

fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    filterLinSpace, horizontalGradientGradient, color='royalblue', label='Gradient-Gradient'
)
ax.plot(
    filterLinSpace, horizontalGradientRandom, color='darkviolet', label='Gradient-Noise'
)
ax.plot(
    filterLinSpace, horizontalGradientLine, color='mediumspringgreen', label='Gradient-Line'
)
plt.setp(ax, xlabel='Arbitrary' , ylabel='Filter Amplitude')
ax.set_title('HU Match Filter')
ax.legend()
ax.grid()

ax = axs[1]
ax.plot(
    filterLinSpace, verticalGradientGradient, color='royalblue', label='Gradient-Gradient'
)
ax.plot(
    filterLinSpace, verticalGradientRandom, color='darkviolet', label='Gradient-Noise'
)
ax.plot(
    filterLinSpace, verticalGradientLine, color='mediumspringgreen', label='Gradient-Line'
)
plt.setp(ax, xlabel='Arbitrary')
ax.set_title('VU Match Filter')
ax.legend()
ax.grid()

#%% Match filtering: 2D Examples 유

# Performing 2D convolution on example maps
convolvedGradient = signal.convolve2d(
    gradientValues,gradientValues,mode='full' #model,data
    )
convolvedRandom = signal.convolve2d(
    gradientValues,randomValues,mode='full'
    )
convolvedLine = signal.convolve2d(
    gradientValues,lineValues,mode='full'
    )

# Plotting results of convolution
fig, axs = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.imshow(convolvedGradient)
ax.set_title('Gradient-Gradient')
ax.grid()

ax = axs[1]
ax.imshow(convolvedRandom)
ax.set_title('Gradient-Noise')
ax.grid()

ax = axs[2]
lastPlot = ax.imshow(convolvedLine)
ax.set_title('Gradient-Line')
ax.grid()

plt.colorbar(lastPlot, fraction=0.050, ax=ax, label = 'δT(z)', pad=0.1)
plt.savefig('2DConvolution.png', dpi=400)

# Plotting unfolded results of convolution to compare for different orientations
horizontalUnfoldedConvolvedGradient = unfolder(convolvedGradient,'horizontal')
horizontalUnfoldedConvolvedRandom = unfolder(convolvedRandom,'horizontal')
horizontalUnfoldedConvolvedLine = unfolder(convolvedLine,'horizontal')
verticalUnfoldedConvolvedGradient = unfolder(convolvedGradient,'vertical')
verticalUnfoldedConvolvedRandom = unfolder(convolvedRandom,'vertical')
verticalUnfoldedConvolvedLine = unfolder(convolvedLine,'vertical')

unfoldedConvolutionSize = len(horizontalUnfoldedConvolvedGradient)
unfoldedConvolutionRange = np.linspace(0,unfoldedConvolutionSize,unfoldedConvolutionSize)

fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    unfoldedConvolutionRange,horizontalUnfoldedConvolvedGradient,color='royalblue',label='Gradient'
)
ax.plot(
    unfoldedConvolutionRange,horizontalUnfoldedConvolvedRandom,color='darkviolet',label='Noise'
)
ax.plot(
    unfoldedConvolutionRange,horizontalUnfoldedConvolvedLine,color='mediumspringgreen',label='Line'
)
ax.set_title('2D Convoltuion HU')
ax.grid()
ax.legend()

ax = axs[1]
ax.plot(
    unfoldedConvolutionRange,verticalUnfoldedConvolvedGradient,color='royalblue',label='Gradient'
)
ax.plot(
    unfoldedConvolutionRange,verticalUnfoldedConvolvedRandom,color='darkviolet',label='Noise'
)
ax.plot(
    unfoldedConvolutionRange,verticalUnfoldedConvolvedLine,color='mediumspringgreen',label='Line'
)
ax.set_title('2D Convoltuion VU')
ax.grid()
ax.legend()


#%% Match filtering: 1D convultion of wake signal 유

# Defining relevant model maps
gradientPositions = wedgeGradient[0] #only within wedge
gradientTemperatures = wedgeGradient[1] #only within wedge
gradientFullPositions = wedgeGradient[2] #the whole space including wedge
gradientFullTemperatures = wedgeGradient[3] #the whole space including wedge
gradientRandom = np.random.random(len(gradientPositions)) #NOT currently Gaussian
gradientFullRandom = np.random.random(len(gradientFullPositions)) #NOT currently Gaussian
scaledGradientPositions = rotatedWedgeGradient[0]
scaledGradientTemperatures = rotatedWedgeGradient[1]
scaledGradientFullPositions = rotatedWedgeGradient[2]
scaledGradientFullTemperatures = rotatedWedgeGradient[3]
scaledGradientRandom = np.random.random(len(scaledGradientPositions)) #NOT currently Gaussian
scaledGradientFullRandom = np.random.random(len(scaledGradientFullPositions)) #NOT currently Gaussian

# Performing 1D match filtering
gradient3DFilter = np.correlate(
    np.ravel(gradientFullTemperatures),np.ravel(gradientFullTemperatures),mode='full'
    )
random3DFilter = np.correlate(
    np.ravel(gradientFullRandom),np.ravel(gradientFullTemperatures),mode='full'
    )

# Plotting results of convolution
fig, axs = plt.subplots(2, 1, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    np.linspace(0,len(gradient3DFilter),len(gradient3DFilter)),gradient3DFilter,
    color='darkviolet',label='Wake-Wake'
)
ax.set_title('Wake Temperature Match Filter')
plt.setp(ax, ylabel='Filter Amplitude')
ax.grid()
ax.legend()

ax = axs[1]
ax.plot(
    np.linspace(0,len(random3DFilter),len(random3DFilter)),random3DFilter,
    color='royalblue',label='Wake-Noise'
)
plt.setp(ax, ylabel='Filter Amplitude', xlabel = 'Arbitrary')
ax.grid()
ax.legend()


#%% Computation of ΛCDM gaussian perturbation noise 유
# Credits to Matteo Blamart

#region
taille=200
fichier1 = open("background_1.txt", "r")
box=np.zeros((taille,taille,taille),float)
li1=[]
k1=0
l=0
for ligne in fichier1:
    l+=1
    a=l//taille
    b=l%taille
    if (l==39999):
        print(li1)
        for z in range(taille):
            li1[z]=float(li1[z])
            box[a][b][z]=li1[z]
        li1 = list(ligne.split(";"))
    else:
        li1 = list(ligne.split(";"))
        del li1[taille]
        for z in range(taille):
            li1[z]=float(li1[z])
            box[a][b][z]=li1[z] #box with first raw=redshift 2nd x axis ans 3rd y axis 
fichier1.close()
#endregion


#%% Algorithm of converting grid points to lattice sheets for line of sight calculations 유
# Everything is currently in Mpc

# Converting position array to position matrix for different z
def arrayToMatrix(array):
    
    # Copying array to no overwrite its values
    tempArray = np.copy(array)
    
    # Finding all unique z values within the array
    tempZValues = []
    for i in range(len(tempArray)):
        tempZValues.append(tempArray[i][-1])
    allZValues = np.array(tempZValues) #all z values within the array
    uniqueZValues = np.unique(allZValues) #unique z values
    numberOfZValues = len(uniqueZValues) #total number of unique z values

    # Defining array points
    newPhysicalPositions = []
    
    # Arranging array elements based on same z value
    for j in uniqueZValues:
        for i in range(len(tempArray)):
            if tempArray[i][-1] == j:
                newPhysicalPositions.append(tempArray[i])
    finalPhysicalPositions = np.array(newPhysicalPositions)  
    
    # Reshaping array to matrix of different z values
    positionsMatrix = finalPhysicalPositions.reshape(numberOfZValues,numberOfZValues**2,3)
      
    return positionsMatrix #matrix with first index for z values

# Defining points for scaled and unscaled wake (will be repeated for scaled wake)
gridPoints = wedgeGradient[2] #(24389, 3)
gridTemps = wedgeGradient[3] #(24389, 1)

# Defining reshaped points
reshapedGridPoints = arrayToMatrix(gridPoints) #(29, 841, 3)
reshapedTemps = gridTemps.reshape(29,841,1) #(29, 841, 1)


#%% Combining wake map with the noise map 유

# ΛCDM gaussian perturbations
perturbationNoise = box*0.001 #converting from mK to K

# Removing a 29x29x29 portion of the noise (scale invariant)
firstDimensionSlice = perturbationNoise[29:58]
secondDimensionSlice = firstDimensionSlice[:,29:58]
thirdDimensionSlice = secondDimensionSlice[:,:,29:58]
# firstDimensionSlice = perturbationNoise[0:29]
# secondDimensionSlice = firstDimensionSlice[:,0:29]
# thirdDimensionSlice = secondDimensionSlice[:,:,0:29]

smallNoiseMap = np.copy(thirdDimensionSlice) #29x29x29 array of temps
reshapedTempMap = gridTemps.reshape((29,29,29)) #29x29x29 array of temps

# Adding both maps
combinedMap = smallNoiseMap + reshapedTempMap

# Showing combination of maps: noise, then wake, then combined
for i in range(29):
    plt.imshow(smallNoiseMap[i])
    plt.show()
    
for i in range(29):
    plt.imshow(reshapedTempMap[i])
    plt.show()

for i in range(29):
    plt.imshow(combinedMap[i])
    plt.show()


#%% Selecting specific slices for plotting purposes 유

# Plotting noise maps
fig, axs = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.imshow(smallNoiseMap[8])
ax.set_title('Slice 9')

ax = axs[1]
ax.imshow(smallNoiseMap[9])
ax.set_title('Slice 10')

ax = axs[2]
lastPlot = ax.imshow(smallNoiseMap[10])
ax.set_title('Slice 11')

plt.colorbar(lastPlot, fraction=0.050, ax=ax, label = 'δT(z)', pad=0.1)
plt.savefig('NoiseMaps.png', dpi=400)


# Plotting wedge maps
fig, axs = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.imshow(reshapedTempMap[8])
ax.set_title('Slice 9')

ax = axs[1]
ax.imshow(reshapedTempMap[9])
ax.set_title('Slice 10')

ax = axs[2]
lastPlot = ax.imshow(reshapedTempMap[10])
ax.set_title('Slice 11')

plt.colorbar(lastPlot, fraction=0.050, ax=ax, label = 'δT(z)', pad=0.1)
plt.savefig('WedgeMaps.png', dpi=400)


# Plotting combined maps
fig, axs = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.imshow(combinedMap[8])
ax.set_title('Slice 9')

ax = axs[1]
ax.imshow(combinedMap[9])
ax.set_title('Slice 10')

ax = axs[2]
lastPlot = ax.imshow(combinedMap[10])
ax.set_title('Slice 11')

plt.colorbar(lastPlot, fraction=0.050, ax=ax, label = 'δT(z)', pad=0.1)
plt.savefig('CombinedMaps.png', dpi=400)

#%% Performing 2D match filtering on the xz-horizontal maps slices for combined maps 유

# Defining example map slices for i = 21 slice
wakeMapModel = reshapedTempMap[15] #0.30 K total
noiseMapData = smallNoiseMap[15] #-1963.74 K total
combinedMapData = combinedMap[15] #-1963.44 K total

# Performing 2D convolution on example maps
convolvedWake = signal.convolve2d(
    wakeMapModel,wakeMapModel,mode='full' #model,data
    )
convolvedNoise = signal.convolve2d(
    wakeMapModel,noiseMapData,mode='full'
    )
convolvedCombined = signal.convolve2d(
    wakeMapModel,combinedMapData,mode='full'
    )

# Plotting results of convolution
fig, axs = plt.subplots(1, 3, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.imshow(convolvedWake)
ax.set_title('Gradient-Gradient')
ax.grid()

ax = axs[1]
ax.imshow(convolvedNoise)
ax.set_title('Gradient-Noise')
ax.grid()

ax = axs[2]
lastPlot = ax.imshow(convolvedCombined)
ax.set_title('Gradient-Combined')
ax.grid()

plt.colorbar(lastPlot, fraction=0.050, ax=ax, label = 'δT(z)', pad=0.1)
plt.savefig('2DConvolution.png', dpi=400)

# Plotting unfolded results of convolution to compare for different orientations
horizontalUnfoldedConvolvedWake = unfolder(convolvedWake,'horizontal')
horizontalUnfoldedConvolvedNoise = unfolder(convolvedNoise,'horizontal')
horizontalUnfoldedConvolvedCombined = unfolder(convolvedCombined,'horizontal')
verticalUnfoldedConvolvedWake = unfolder(convolvedWake,'vertical')
verticalUnfoldedConvolvedNoise = unfolder(convolvedNoise,'vertical')
verticalUnfoldedConvolvedCombined = unfolder(convolvedCombined,'vertical')

unfoldedConvolutionSize = len(horizontalUnfoldedConvolvedWake)
unfoldedConvolutionRange = np.linspace(0,unfoldedConvolutionSize,unfoldedConvolutionSize)

fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    unfoldedConvolutionRange,horizontalUnfoldedConvolvedWake,color='royalblue',label='Gradient'
)
ax.plot(
    unfoldedConvolutionRange,horizontalUnfoldedConvolvedNoise,color='darkviolet',label='Noise'
)
ax.plot(
    unfoldedConvolutionRange,horizontalUnfoldedConvolvedCombined,color='mediumspringgreen',label='Combined'
)
ax.set_title('2D Convoltuion HU')
ax.grid()
ax.legend()

ax = axs[1]
ax.plot(
    unfoldedConvolutionRange,verticalUnfoldedConvolvedWake,color='royalblue',label='Gradient'
)
ax.plot(
    unfoldedConvolutionRange,verticalUnfoldedConvolvedNoise,color='darkviolet',label='Noise'
)
ax.plot(
    unfoldedConvolutionRange,verticalUnfoldedConvolvedCombined,color='mediumspringgreen',label='Combined'
)
ax.set_title('2D Convolution VU')
ax.grid()
ax.legend()


#%% Performing 1D match filtering on the xz-horizontal maps slices for combined maps 유

# Unfolding maps for match filtering for i = 21 slice
horizontalUnfoldedNoise = unfolder(noiseMapData,'horizontal')
verticalUnfoldedNoise = unfolder(noiseMapData,'vertical')
horizontalUnfoldedWake = unfolder(wakeMapModel, 'horizontal')
verticalUnfoldedWake = unfolder(wakeMapModel, 'vertical')
horizontalUnfoldedCombined = unfolder(combinedMapData,'horizontal')
verticalUnfoldedCombined = unfolder(combinedMapData,'vertical')

unfoldSize = len(horizontalUnfoldedWake)
unfoldedIndices = np.linspace(0,unfoldSize,unfoldSize)
indicesSize = len(unfoldedIndices)/2

# Plotting results of unfolding (without noise)
fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    unfoldedIndices, horizontalUnfoldedWake, color='royalblue', label='Gradient'
)
ax.plot(
    unfoldedIndices, horizontalUnfoldedCombined, color='mediumspringgreen', label='Combined'
)
plt.setp(ax, xlabel='Index' , ylabel='Brightness Temperature [K]')
ax.set_ylim(-0.015,0.015)
ax.set_xlim(indicesSize-200,indicesSize+200)
ax.set_title('Horizontal Unfolding')
ax.legend()
ax.grid()

ax = axs[1]
ax.plot(
    unfoldedIndices, verticalUnfoldedWake, color='royalblue', label='Gradient'
)
ax.plot(
    unfoldedIndices, verticalUnfoldedCombined, color='mediumspringgreen', label='Combined'
)
plt.setp(ax, xlabel='Index')
ax.set_ylim(-0.015,0.015)
ax.set_xlim(indicesSize-200,indicesSize+200)
ax.set_title('Vertical Unfolding')
ax.legend()
ax.grid()

# Plotting results of unfolding (with noise)
fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    unfoldedIndices, horizontalUnfoldedNoise, color='darkviolet', label='Noise'
)
plt.setp(ax, xlabel='Index' , ylabel='Brightness Temperature [K]')
ax.set_xlim(indicesSize-200,indicesSize+200)
ax.set_ylim(-0.0185,0)
ax.set_title('Horizontal Unfolding')
ax.legend()
ax.grid()

ax = axs[1]
ax.plot(
    unfoldedIndices, verticalUnfoldedNoise, color='darkviolet', label='Noise'
)
plt.setp(ax, xlabel='Index')
ax.set_xlim(indicesSize-200,indicesSize+200)
ax.set_ylim(-0.0185,0)
ax.set_title('Vertical Unfolding')
ax.legend()
ax.grid()

# Performing match filtering on the data sample
horizontalWakeWake = np.correlate(
    horizontalUnfoldedConvolvedWake,horizontalUnfoldedConvolvedWake,mode='full'
)
horizontalWakeNoise = np.correlate(
    horizontalUnfoldedConvolvedWake,horizontalUnfoldedConvolvedNoise,mode='full'
)
horizontalWakeCombined = np.correlate(
    horizontalUnfoldedConvolvedWake,horizontalUnfoldedConvolvedCombined,mode='full'
)
verticalWakeWake = np.correlate(
    verticalUnfoldedConvolvedWake,verticalUnfoldedConvolvedWake,mode='full'
)
verticalWakeNoise = np.correlate(
    verticalUnfoldedConvolvedWake,verticalUnfoldedConvolvedNoise,mode='full'
)
verticalWakeCombined = np.correlate(
    verticalUnfoldedConvolvedWake,verticalUnfoldedConvolvedCombined,mode='full'
)

matchFilterSize = len(horizontalWakeWake)
filterLinSpace = np.linspace(0,matchFilterSize,matchFilterSize)

# Plotting results of 1D match filtering
fig, axs = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

ax = axs[0]
ax.plot(
    filterLinSpace, horizontalWakeWake, color='royalblue', label='Gradient-Gradient'
)
ax.plot(
    filterLinSpace, horizontalWakeNoise, color='darkviolet', label='Gradient-Noise'
)
ax.plot(
    filterLinSpace, horizontalWakeCombined, color='mediumspringgreen', label='Gradient-Combined'
)
plt.setp(ax, xlabel='Arbitrary' , ylabel='Filter Amplitude')
ax.set_title('HU Match Filter')
ax.legend()
ax.grid()

ax = axs[1]
ax.plot(
    filterLinSpace, verticalWakeWake, color='royalblue', label='Gradient-Gradient'
)
ax.plot(
    filterLinSpace, verticalWakeNoise, color='darkviolet', label='Gradient-Noise'
)
ax.plot(
    filterLinSpace, verticalWakeCombined, color='mediumspringgreen', label='Gradient-Combined'
)
plt.setp(ax, xlabel='Arbitrary')
ax.set_title('VU Match Filter')
ax.legend()
ax.grid()


#%% CODE GRAVEYARD 유

#region

###Testing created convolution functions to build up to 3D 유

# Defining 1D convolution: s_t = h_(t-k) * d^k, summed over k from -n^2/2 to n^2/2 
# def matchFilter(data,model):
#     size = len(data)
#     # sumLimits = (size**2)/2
#     s = np.zeros(size)
#     for t in range(size): #for different t values
#         for k in range(-int(size),int(size)): #-n^2/2 to n^2/2
#             s[t] += model[(t-k)%size]*data[(k)%size]
#     return s

# # Testing defined function on example maps & plotting
# hGradientGradient = matchFilter(horizontalGradientGradient,horizontalGradientGradient)
# hGradientRandom = matchFilter(horizontalGradientGradient,horizontalGradientRandom)
# hGradientLine = matchFilter(horizontalGradientGradient,horizontalGradientLine)
# mFilterSize = len(hGradientGradient)
# fLinSpace = np.linspace(0,mFilterSize,mFilterSize)

# plt.plot(fLinSpace,hGradientGradient,color='royalblue', label='Gradient-Gradient')
# plt.plot(fLinSpace,hGradientRandom,color='darkviolet', label='Gradient-Noise')
# plt.plot(fLinSpace,hGradientLine,color='mediumspringgreen', label='Gradient-Line')
# plt.title('HU Manual Match Filter')
# plt.xlabel('Arbitrary')
# plt.ylabel('Filter Amplitude')
# plt.legend()
# plt.grid()


### Algorithm of converting grid points to lattice sheets for line of sight calculations 유

# # Creating n x n mesh grid
# n = 3
# meshGrid = np.mgrid[0:n:1,0:n:1,0:n:1]
# X, Y, Z = meshGrid
# positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])

# # Converting mesh grid to (n^2,3) array of physical points
# emptyPositions = []
# for i in range(len(positions.T)):
#     emptyPositions.append([positions[0][i],positions[1][i],positions[2][i]])
# physicalPositions = np.array(emptyPositions)

# # Converting position array to position matrix for different z
# def arrayToMatrix(array):
    
#     # Defining array points
#     tempArray = np.copy(array)
#     newPhysicalPositions = []
#     n = int(np.cbrt(len(array))) #getting n for an nxnxn matrix from an (n^3,3) matrix
    
#     # Arranging array elements based on same z value
#     for j in range(n):
#         for i in range(len(tempArray)):
#             if tempArray[i][-1] == j:
#                 newPhysicalPositions.append(tempArray[i])
#     finalPhysicalPositions = np.array(newPhysicalPositions)  
    
#     # Reshaping array to matrix of different z values
#     positionsMatrix = finalPhysicalPositions.reshape(n,n**2,3)
      
#     return positionsMatrix #matrix with first index for z values

# # # Example converted matrix for n x n
# # latticePoints = arrayToMatrix(physicalPositions)
# # latticeTemperatures = np.arange(0,n*(n**2),1).reshape((n,n**2))
# # print(latticePoints)
# # print(latticeTemperatures)

# # # Now creating function to convert Matteo's 200x200x200 lattice to a matrix to work with
# # testLattice = np.array([1,2,3,4,5,6,7,8]).reshape((2,2,2))

#endregion

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

> [!TIP]
> For a detailed derivation of the brightness temperature formula, the homotopy classification of defects, and scaling solutions for cosmic strings, consult the PDF in the main repository.

>[!NOTE]
> Once the method is validated on simulated data, it can be applied to actual 21cm observations (e.g., from upcoming radio interferometers) to search for cosmic string imprints in the real universe.
