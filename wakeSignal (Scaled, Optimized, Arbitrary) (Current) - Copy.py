#%% Importing modules & notes 유
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from astropy.cosmology import Planck18
from astropy.cosmology import z_at_value
import astropy.units as u
import random as rd


#%% Defining compact plotting functions 유

# Function that quickly plots 3D figures
def quickPlot(scale,points,labels,colors,colormaps=[],redshiftScale=False):
    
    # Defining some plotting parameters and fixed scales 
    fig = plt.figure(figsize=(12, 10), dpi=80)
    ax = fig.add_subplot(111, projection="3d")
    redshiftScale = 6.77E-5
    
    # Finding the approximate center point of wake to scale the plots properly
    approximateCenter = []
    for i in points:
        if len(i) == 6:
            for j in range(0,3):
                approximateCenter.append(
                    0.5*(np.max(i.T[j])-np.min(i.T[j]))
                )

    # Labelling & scaling axes depending if physical space is used, or 2D x redshift space
    if redshiftScale == False: 
        
        ax.set_xlabel('x [Mpc]')
        ax.set_ylabel('y [Mpc]')
        ax.set_zlabel('z [Mpc]')
        ax.set_xlim(approximateCenter[0]-scale/2,approximateCenter[0]+scale/2)
        ax.set_ylim(approximateCenter[1]-scale/2,approximateCenter[1]+scale/2)
        ax.set_zlim(approximateCenter[2]-scale/2,approximateCenter[2]+scale/2)
    
    if redshiftScale == True:
        
        ax.set_xlabel('x [Mpc]')
        ax.set_ylabel('y [Mpc]')
        ax.set_zlabel('z')
        ax.set_xlim(approximateCenter[0]-scale/2,approximateCenter[0]+scale/2)
        ax.set_ylim(approximateCenter[1]-scale/2,approximateCenter[1]+scale/2)
        ax.set_zlim(approximateCenter[2]-redshiftScale/2,approximateCenter[2]+redshiftScale/2)
    
    # Differentiating if points have an optional temperature color map
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
    else:
        physicalPlot = ax.scatter(
            points[:,0],points[:,1],points[:,2], c = colormaps,
            label = 'Internal temperature points', cmap = plt.cm.viridis, s = 90
        )
    
        fig.colorbar(physicalPlot, ax=ax, label = 'δT(z)', pad=0.1)
    
    ax.legend()
    plt.show()

    return


# Function that quickly plots 3D figures with convex hulls
def quickPlotHull(scale,hullVertices, testPoints=[],redshiftScale = False):
    
    # Defining some plotting parameters and fixed scales 
    fig = plt.figure(figsize=(12, 10), dpi=80)
    ax = fig.add_subplot(111, projection="3d")
    redshiftScale = 6.77E-5
    
    # Finding the approximate center point of wake to scale the plots properly
    approximateCenter = []
    for i in hullVertices:
        if len(i) == 6:
            for j in range(0,3):
                approximateCenter.append(
                    0.5*(np.max(i.T[j])-np.min(i.T[j]))
                )

    # Labelling & scaling axes depending if physical space is used, or 2D x redshift space
    if redshiftScale == False: 
        
        ax.set_xlabel('x [Mpc]')
        ax.set_ylabel('y [Mpc]')
        ax.set_zlabel('z [Mpc]')
        ax.set_xlim(approximateCenter[0]-scale/2,approximateCenter[0]+scale/2)
        ax.set_ylim(approximateCenter[1]-scale/2,approximateCenter[1]+scale/2)
        ax.set_zlim(approximateCenter[2]-scale/2,approximateCenter[2]+scale/2)
    
    if redshiftScale == True:
        
        ax.set_xlabel('x [Mpc]')
        ax.set_ylabel('y [Mpc]')
        ax.set_zlabel('z')
        ax.set_xlim(approximateCenter[0]-scale/2,approximateCenter[0]+scale/2)
        ax.set_ylim(approximateCenter[1]-scale/2,approximateCenter[1]+scale/2)
        ax.set_zlim(approximateCenter[2]-redshiftScale/2,approximateCenter[2]+redshiftScale/2)
    
    # Defining hull & plotting hull vertices
    hull = ConvexHull(hullVertices)
    ax.scatter(
        hullVertices.T[0], hullVertices.T[1], hullVertices.T[2], color='red', s = 90
        )
    
    # Plotting hull simplices that are connected to vertices
    for s in hull.simplices:
        s = np.append(s, s[0])
        ax.plot(hullVertices[s, 0], hullVertices[s, 1], hullVertices[s, 2], color='maroon')
    
    # Plotting any potential test points along side the hull
    if len(testPoints) != 0:
        ax.scatter(
        testPoints[:,0],testPoints[:,1],testPoints[:,2], color='black', s = 90
        )
    
    plt.show()

    return


# %% Defining Hubble volume and cosmic string wake wedge 유

# Initializing some required parameters
hubbleParameter = 70.41*1000 #[m/s/Mpc] H(z=0)=H0
densityRatio = 0.3 #non-relativistic matter density for flat ΛCDM 
speedOfLight = 299792458 #[m/s]
hubbleScale = speedOfLight/hubbleParameter #[Mpc] Hubble length
hubbleLattice = hubbleScale * np.array([
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
     ])

# Defining various required parameters
Gμ = 1E-8 #string tension
deficitAngle = np.pi/3 #actually is 8*np.pi*Gμ, using π/3 to visualize angle
c1 = rd.uniform(1-0.01,1+0.01) #O(1)
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
    ['blue','red']
    )

# Plotting again with smaller axis
quickPlot(
    hubbleScale/scaleConverter,[wakeWedge, wakeMiddlePoint],['Wake wedge','Wake center'],
    ['red','black']
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
    ['darkorchid','red','black']
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
ax.legend()
ax.grid()

# Redshift as a function of distance plot
ax = axs[1]
ax.plot(y, distanceToRedshift(y), color='red', label='z(d)')
ax.grid()
ax.legend()

# Demonstrating inverse function equivalence
print('d(z=100): ', int(redshiftToDistance(100)), 'Mpc')
print('z(d~12831 Mpc): ', int(distanceToRedshift(12831.379081434678)))

# Defining z(d) scaling function to apply to arbitrary wake wedges
def wedgeToRedshift(wedge):
    
    wedgeCopy = np.copy(wedge)
    
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
quickPlot(
    hubbleScale/scaleConverter,[scaledRotatedWedge, scaledWedge, scaledWedgeCenter],['Rotated wake','Wake wedge','Wake center'],
    ['darkorchid','red','black'],redshiftScale=True
    )


#%% Computations to convert comoving coordinates to expanding physical coordinates

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
quickPlot(
    hubbleScale/scaleConverter,[physicalScaledRotatedWedge, physicalScaledWedge, physicalScaledWedgeCenter],
    ['Physical rotated wake','Physical Wake wedge','Physical wake center'],
    ['darkorchid','red','black'],redshiftScale=True
    )


#%% Computing the boundary of the wedge via convex hulls & calculating if test points are within it 유

# Using cube as an initial example (unscaled Hubble lattice)
cubeVertices = hubbleScale/scaleConverter* np.array([
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
     ])

# Plotting cube complex closure
quickPlotHull(hubbleScale/scaleConverter,cubeVertices)

# Repeating the above calculation with the wake wedge
quickPlotHull(hubbleScale/scaleConverter,wakeWedge)
quickPlotHull(hubbleScale/scaleConverter,physicalScaledWedge, redshiftScale=True)
quickPlotHull(hubbleScale/scaleConverter,physicalScaledRotatedWedge, redshiftScale=True)

# Function that checks if a set of points lie within the hull
def pointChecker(scale,hullPoints, testPoints, plot=False,redshiftScale=False):
    
    # Defining hull and test points
    hull = ConvexHull(hullPoints)
    newPoints = np.append(hullPoints, testPoints, axis=0)
    newHull = ConvexHull(newPoints)
    
    # Checking if points are within the hull
    if list(hull.vertices) == list(newHull.vertices):
        criterion = True
    else:
        criterion = False
    
    # Whether or not this function shows the plots
    if plot == True:
    
        # Plotting the hull & test points
        quickPlotHull(scale,hullPoints,testPoints,redshiftScale)
    
    else:
        
        return criterion
    
# Running preliminary examples
pointChecker(hubbleScale/scaleConverter,wakeWedge,np.array([[0.1,0.15,0.05]]))
print(pointChecker(hubbleScale/scaleConverter,wakeWedge,np.array([[0.1,0.15,0.05]])))

pointChecker(hubbleScale/scaleConverter,physicalScaledRotatedWedge,np.array([[0.2,0.1,1E-5]]),redshiftScale=True)
print(pointChecker(hubbleScale/scaleConverter,physicalScaledRotatedWedge,np.array([[0.2,0.1,1E-5]]),redshiftScale=True))


#%% Computing brightness temperature function δT(z) for points within wake 유

# Defining functions for brightness temperature function
kineticTemperature = lambda z : (22020*Gμ*(1E12)*(stringSpeed*gammaFactor)**2*speedOfLight**(-2))/(z+1) #[K]
deexcitationCrossSection = lambda z : 2.76E-9 * (1 - 1/((1+0.0007112 * (kineticTemperature(z))**(2.28) )**(0.014)))*(1/(100**3))
photonCMBTemp = lambda z : 2.72548*(1+z) #[K] from wiki
densityNumber = lambda z : 0.76* (1+z)**3 #[1/m^3]
collisionCoefficient = lambda z : (densityNumber(z)*deexcitationCrossSection(z)*0.068)/(2.80E-15 * photonCMBTemp(z))
brightnessTemperature = lambda z : 0.07 * (collisionCoefficient(z)/(1+collisionCoefficient(z))) * (1-photonCMBTemp(z)/kineticTemperature(z)) * np.sqrt(1+z) #[K] δT(z)
brightnessDistanceTemperature = lambda d : brightnessTemperature(distanceToRedshift(d)) #[K] δT(d)

# Defining function to find points within wedge and apply temperature gradient
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
            approximateCenter[2]-scale/2:approximateCenter[2]+scale/2:resolution #was 0.3 and 0.01
            ] #in each axis it ranges from 0 to 0.3 in 0.01 steps
        X, Y, Z = meshGrid
        positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
    
    else:
        
        # Finding grid of points within the Hubble lattice
        meshGrid = np.mgrid[
            approximateCenter[0]-scale/2:approximateCenter[0]+scale/2:resolution, 
            approximateCenter[1]-scale/2:approximateCenter[1]+scale/2:resolution, 
            approximateCenter[2]-redshiftScale/2:approximateCenter[2]+redshiftScale/2:1.67E-6 #was 0.3 and 0.01
            ] #in each axis it ranges from 0 to 0.3 in 0.01 steps
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
    
    return (physicalWedgePoints,physicalWedgeTemperatures)

wedgeGradient = gradientManifold(hubbleScale/scaleConverter,0.01,wakeWedge)
rotatedWedgeGradient = gradientManifold(hubbleScale/scaleConverter,0.01,physicalScaledRotatedWedge,redshiftMpc=[9])

# Plotting physical points within wake excluding ambient points
quickPlot(
    hubbleScale/scaleConverter,wedgeGradient[0],['Internal temperature points'],
    ['red'], wedgeGradient[1]
    )

quickPlot(
    hubbleScale/scaleConverter,rotatedWedgeGradient[0],['Internal temperature points'],
    ['red'], rotatedWedgeGradient[1], redshiftScale=True
    )


#%% Defining ΛCDM gaussian perturbations 유


#%% Superimposing wake signal and perturbations in a plot 유


#%% Computing 3D match filtering statistic 유