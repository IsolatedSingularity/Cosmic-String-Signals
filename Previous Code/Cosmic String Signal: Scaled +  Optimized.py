#%% Importing modules & notes 유
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from astropy.cosmology import Planck18
from astropy.cosmology import z_at_value
import astropy.units as u


#%% Task list 유 
# Testing convex hull: plot temperature points and hull
# Make fixed values parameters
# Wake length dimensions
# Scale dimensions from Mpc to redshift (non-natural, astrophysical units) (take amplitude in 21cmfast and convert to natural untis)
# Scale d -> z(d) from string theory
# Combine: Resolution, frequency (21cmFAST Mpc units, output amplitude in Kelvin) & match grid points (EXPLICIT), Gmu free param (start with 10^-8)
# Clean up all code and plots & optimize
# Analytically how does crossover redshift depend on Gmu
#z~recombination
#different resolutions in z and angular directions (angular resolution can be low)
#width of wake is growing in comoving coords (double split figure) (later t thicker wake)
#time interval
#match filtering -> other statistics (depending on Gmu)


#%% Defining compact plotting functions 유

# Function that quickly plots 3D figures
def quickPlot(scale,points,labels,colors,colormaps=[],redshiftScale = False):
    
    # Defining generic plotting parameters 
    fig = plt.figure(figsize=(12, 10), dpi=80)
    ax = fig.add_subplot(111, projection="3d")

    # Labelling & scaling axes depending if physical space is used, or 2D x redshift space
    if redshiftScale == False: 
        
        for axis in ["x", "y", "z"]:
            scale = scale #initializing local scaling variable to be used by eval()
            eval("ax.set_{:s}label('{:s} [Mpc]')".format(axis, axis)) #labelling axes by x,y,z
            eval("ax.set_{:s}lim(0,scale)".format(axis)) #setting axes limits as (0,scale)
    
    else:
        
        for axis in ["x", "y"]:
            scale = scale #initializing local scaling variable to be used by eval()
            eval("ax.set_{:s}label('{:s} [Mpc]')".format(axis, axis)) #labelling axes by x,y,z
            eval("ax.set_{:s}lim(0,scale)".format(axis)) #setting axes limits as (0,scale)
        ax.set_zlabel('z')
        ax.set_zlim(0,4.5E-5)
    
    # Differentiating if points have an optional color map
    if len(colormaps) == 0:
        
        # Iterating through sets of points
        for pointType in range(len(points)):
            
            # Differentiating between individual and sets of points
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
    
        # ax.set_box_aspect([np.ptp(i) for i in meshGrid]) #equal aspect ratio for axes
        fig.colorbar(physicalPlot, ax=ax, label = 'δT(z)', pad=0.1)
    
    ax.legend()
    plt.show()

    return


# Labelling & scaling axes depending if physical space is used, or 2D x redshift space
def quickPlotHull(scale,hullVertices, testPoints=[],redshiftScale = False):
    
    # Defining generic plotting parameters
    fig = plt.figure(figsize=(12, 10), dpi=80)
    ax = fig.add_subplot(111, projection="3d")
    
    # Labelling axes
    if redshiftScale == False: 
        
        for axis in ["x", "y", "z"]:
            scale = scale #initializing local scaling variable to be used by eval()
            eval("ax.set_{:s}label('{:s} [Mpc]')".format(axis, axis)) #labelling axes by x,y,z
            eval("ax.set_{:s}lim(0,scale)".format(axis)) #setting axes limits as (0,scale)
    
    else:
        
        for axis in ["x", "y"]:
            scale = scale #initializing local scaling variable to be used by eval()
            eval("ax.set_{:s}label('{:s} [Mpc]')".format(axis, axis)) #labelling axes by x,y,z
            eval("ax.set_{:s}lim(0,scale)".format(axis)) #setting axes limits as (0,scale)
        ax.set_zlabel('z')
        ax.set_zlim(0,4.5E-5)
    
    # Defining hull & plotting hull vertices
    hull = ConvexHull(hullVertices)
    ax.scatter(
        hullVertices.T[0], hullVertices.T[1], hullVertices.T[2], color='red', s = 90
        )
    
    # Plotting hull simplices that are connected to vertices
    for s in hull.simplices:
        s = np.append(s, s[0])  #permuting back to first coordinate
        ax.plot(hullVertices[s, 0], hullVertices[s, 1], hullVertices[s, 2], color='maroon')
    
    # Plotting any potential test points
    if len(testPoints) != 0:
        ax.scatter(
        testPoints[:,0],testPoints[:,1],testPoints[:,2], color='black', s = 90
        )
    
    plt.show()

    return


# %% Defining Hubble volume and cosmic string wake wedge 유

# Defining ambient grid and Hubble lattice/volume
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
c1 = 0.88 #O(1)
meterToMpc = 3.24078E-23 #1m = 3.24078E-23 Mpc; to convert m -> Mpc or rather m/Mpc -> 1 
redshiftTime = lambda z : 2/(3 * hubbleParameter * (densityRatio**(1/2)) * (1+z)**(3/2)) * meterToMpc**(-1) #t(z) in s
formationTime = redshiftTime(1100) #time during which string formation is most prominent: recombination z~1100

# Defining the wake geometry in 3D physical space
stringSpeed = speedOfLight * 0.8 #[m/s] conservative estimate
gammaFactor = ( 1/( 1 - (stringSpeed/speedOfLight)**2 ) )**(1/2)
wakeLength = c1 * formationTime * speedOfLight * meterToMpc #[Mpc]
wakeDepth = formationTime * gammaFactor * stringSpeed * meterToMpc #[Mpc] radial length
wakeTipPoint = np.array([0.06,0.06,0.06])
wakeEndPoints = np.array([
    [wakeTipPoint[0]+wakeDepth*np.cos(deficitAngle/2),wakeTipPoint[1]+wakeDepth*np.sin(deficitAngle/2),wakeTipPoint[2]],
    [wakeTipPoint[0]+wakeDepth*np.cos(deficitAngle/2),wakeTipPoint[1]-wakeDepth*np.sin(deficitAngle/2),wakeTipPoint[2]]
    ])
projectedWakePoints = np.array([wakeTipPoint+[0,0,wakeLength],wakeEndPoints[0]+[0,0,wakeLength],wakeEndPoints[1]+[0,0,wakeLength]])
middlePointWake = np.array([wakeTipPoint[0]+(1/2)*wakeDepth*np.cos(deficitAngle/2),wakeTipPoint[1],wakeTipPoint[2]+wakeLength/2])
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
    hubbleScale/scaleConverter,[hubbleLattice, wakeWedge, middlePointWake],['Hubble lattice','Wake wedge','Wake center'],
    ['blue','red','black']
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
shiftedWedge = np.add(wakeWedge,-middlePointWake)

# Applying rotations on wake (random angles around all axes)
modifiedWedgeI = np.dot(shiftedWedge, xAxisRotation(np.pi/2).T) #rotation of pi/2 about x-axis
modifiedWedgeII = np.dot(modifiedWedgeI, yAxisRotation(np.pi/3).T) #rotation of pi/3 about y-axis
modifiedWedgeIII = np.dot(modifiedWedgeII, zAxisRotation(2*np.pi).T) #rotation of 2 pi about z-axis

# Shifting wake back to center created
rotatedWedge = np.add(modifiedWedgeIII,middlePointWake)

# Plotting unrotated and rotated wake in the Hubble lattice
quickPlot(
    hubbleScale/scaleConverter,[rotatedWedge, wakeWedge, middlePointWake],['Rotated wake','Wake wedge','Wake center'],
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
    if len(wedge) == 3:
        wedgeCopy[-1] = distanceToRedshift(wedge[-1])
        
    else:
        
        for point in range(len(wedge)):
            wedgeCopy[point][-1] = distanceToRedshift(wedge[point][-1])
            
    return wedgeCopy

scaledWedge = wedgeToRedshift(wakeWedge)
scaledRotatedWedge = wedgeToRedshift(rotatedWedge)
scaledWedgeCenter = wedgeToRedshift(middlePointWake)

# Plotting scaled (d [Mpc] -> z) wedges in 2D physical x 1D redshift space aka (2,1) space
quickPlot(
    hubbleScale/scaleConverter,[scaledRotatedWedge, scaledWedge, scaledWedgeCenter],['Rotated wake','Wake wedge','Wake center'],
    ['darkorchid','red','black'],redshiftScale=True
    )


#%% Computations to convert comoving coordinates to expanxing physical coordinates


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
quickPlotHull(hubbleScale/scaleConverter,scaledWedge, redshiftScale=True)
quickPlotHull(hubbleScale/scaleConverter,scaledRotatedWedge, redshiftScale=True)

# Checking if a set of points lie within the hull
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
pointChecker(hubbleScale/scaleConverter,wakeWedge,np.array([[0.1,0.15,0.05]]),True)
print(pointChecker(hubbleScale/scaleConverter,wakeWedge,np.array([[0.1,0.15,0.05]])))

pointChecker(hubbleScale/scaleConverter,scaledRotatedWedge,np.array([[0.2,0.1,1E-5]]),True,redshiftScale=True)
print(pointChecker(hubbleScale/scaleConverter,wakeWedge,np.array([[0.2,0.1,1E-5]]),redshiftScale=True))


#%% Computing brightness temperature function δT(z) for points within wake 유

# Defining temporary function for brightness temperature
distanceBrightnessTemperature = lambda d : np.sqrt(1+distanceToRedshift(d))
redshiftBrightnessTemperature = lambda z : np.sqrt(1+z)

# Defining function to find points within wedge and apply temperature gradient
def gradientManifold(physicalScale,resolution,wedge,redshiftMpc=[]):
    
    # Determining if z axis is scaled to redshift, and if so what physicalScale should be used for the xy plane
    if len(redshiftMpc) == 0:

        # Finding grid of points within the Hubble lattice
        meshGrid = np.mgrid[
            0:physicalScale:resolution, 0:physicalScale:resolution, 0:physicalScale:resolution #was 0.3 and 0.01
            ] #in each axis it ranges from 0 to 0.3 in 0.01 steps
        X, Y, Z = meshGrid
        positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
    
    else:
        
        # Finding grid of points within the Hubble lattice
        minZ = min([row[-1] for row in wedge])
    
        meshGrid = np.mgrid[
            0:redshiftMpc[0]:redshiftMpc[1], 0:redshiftMpc[0]:redshiftMpc[1], minZ:physicalScale:resolution #was 0.3 and 0.01
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
    
    if len(redshiftMpc) == 0:
        for i in range(len(physicalPositions)):
            if pointChecker(hubbleScale/scaleConverter,wakeCopy,np.array([physicalPositions[i]])) == True:
                temperatureValues[i] = distanceBrightnessTemperature(physicalPositions[i][2])
                wedgePoints.append(physicalPositions[i])
                wedgeTemperatures.append(distanceBrightnessTemperature(physicalPositions[i][2]))    
     
    # Determining if z axis is scaled to redshift           
    else:
        for i in range(len(physicalPositions)):
            if pointChecker(hubbleScale/scaleConverter,wakeCopy,np.array([physicalPositions[i]])) == True:
                temperatureValues[i] = redshiftBrightnessTemperature(physicalPositions[i][2])
                wedgePoints.append(physicalPositions[i])
                wedgeTemperatures.append(redshiftBrightnessTemperature(physicalPositions[i][2]))   
        
    physicalWedgePoints = np.array(wedgePoints)
    physicalWedgeTemperatures = np.array(wedgeTemperatures)
    
    return (physicalWedgePoints,physicalWedgeTemperatures)

wedgeGradient = gradientManifold(0.3,0.01,wakeWedge)
rotatedWedgeGradient = gradientManifold(5E-5,1.67E-6,scaledRotatedWedge,redshiftMpc=[0.3,0.01])

# Plotting physical points within wake excluding ambient points
quickPlot(
    hubbleScale/scaleConverter,wedgeGradient[0],['Internal temperature points'],
    ['red'], wedgeGradient[1]
    )

quickPlot(
    hubbleScale/scaleConverter,rotatedWedgeGradient[0],['Internal temperature points'],
    ['red'], rotatedWedgeGradient[1], redshiftScale=True
    )


#%% Combining code with Matteo 유
#region
# #Computing brightness temperature function δT(z) for points within wake 유

# # Defining temporary function for brightness temperature
# def brightnessTemperature(d): #function of position currently
#     return np.sqrt(1+distanceToRedshift(d))

# # Finding grid of points within (unscaled) Hubble lattice (3 x 1000 elements)
# meshGrid = np.mgrid[0:1:0.02, 0:1:0.02, 0:1:0.02] #in each axis it ranges from 0 to 0.3 in 0.01 steps
# X, Y, Z = meshGrid
# positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])

# # Getting array of all points within the space
# emptyPositions = []
# wedgePoints = [] #Physical points within wedge only
# wedgeTemperatures = [] #Temperature of points within wedge only

# # Converting mesh grid to physical points in 3D space
# for i in range(len(positions.T)):
#     emptyPositions.append([positions[0][i],positions[1][i],positions[2][i]])

# physicalPositions = np.array(emptyPositions)

# # Iterating through the Hubble lattice to find points within wake and assign a temperature
# temperatureValues = np.zeros((len(positions.T),1))

# for i in range(len(physicalPositions)):
#     if pointChecker(wakeWedge,np.array([physicalPositions[i]])) == True:
#         temperatureValues[i] = brightnessTemperature(physicalPositions[i][2])
#         wedgePoints.append(physicalPositions[i])
#         wedgeTemperatures.append(brightnessTemperature(physicalPositions[i][2]))
        
# physicalWedgePoints = np.array(wedgePoints)
# physicalWedgeTemperatures = np.array(wedgeTemperatures)

# # Plotting physical points within wake excluding ambient points
# quickPlot(
#     hubbleScale/scaleConverter,physicalWedgePoints,['Internal temperature points'],
#     ['red'], physicalWedgeTemperatures
#     )
#endregion
