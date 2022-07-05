#%% Importing modules & notes 유
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

# Task list
# Testing convex hull: plot temperature points and hull
# Make fixed values parameters
# Wake length dimensions
# Scale dimensions from Mpc to redshift (non-natural, astrophysical units) (take amplitude in 21cmfast and convert to natural untis)
# Scale d -> z(d) from string theory
# Combine: Resolution, frequency (21cmFAST Mpc units, output amplitude in Kelvin) & match grid points (EXPLICIT), Gmu free param (start with 10^-8)
# Clean up all code and plots & optimize
# Analytically how does crossover redshift depend on Gmu


#%% Defining compact plotting functions 유

def quickPlot(scale,points,labels,colors,colormaps=[]):
    
    # Defining generic plotting parameters 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Labelling & scaling axes
    for axis in ["x", "y", "z"]:
        scale = scale #Initializing local scaling variable to be used by eval()
        eval("ax.set_{:s}label('{:s}')".format(axis, axis)) #Labelling axes by x,y,z
        eval("ax.set_{:s}lim(0,scale)".format(axis)) #Setting axes limits as (0,scale)
    
    # Differentiating if points have optional (temperature) color maps
    if len(colormaps) == 0:
        
        # Iterating through sets of points
        for pointType in range(len(points)):
            
            # Differentiating between individual and sets of points
            if points[pointType].ndim == 2:
                ax.scatter(
                    points[pointType][:,0], points[pointType][:,1], points[pointType][:,2],
                    label = labels[pointType], color = colors[pointType]
                )
                
            else:
                ax.scatter(
                    points[pointType][0], points[pointType][1], points[pointType][2],
                    label = labels[pointType], color = colors[pointType]
                )
    else:
        physicalPlot = ax.scatter(
        points[:,0],points[:,1],points[:,2], c = colormaps,
        label = 'Internal temperature points', cmap = plt.cm.viridis
        )
    
        ax.set_box_aspect([np.ptp(i) for i in meshGrid]) #Equal aspect ratio for axes
        fig.colorbar(physicalPlot, ax=ax, label = 'δT(z)', pad=0.1)
    
    ax.legend()
    plt.show()

    return


def quickPlotHull(hullVertices, testPoints=[]):
    
    # Defining generic plotting parameters
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    
    # Labelling axes
    for axis in ["x", "y", "z"]:
        eval("ax.set_{:s}label('{:s}')".format(axis, axis)) #Labelling axes by x,y,z
    
    # Defining hull & plotting hull vertices
    hull = ConvexHull(hullVertices)
    ax.plot(
        hullVertices.T[0], hullVertices.T[1], hullVertices.T[2], "ko", color='blue'
        )
    
    # Plotting hull simplices that are connected to vertices
    for s in hull.simplices:
        s = np.append(s, s[0])  # Permuting back to first coordinate
        ax.plot(hullVertices[s, 0], hullVertices[s, 1], hullVertices[s, 2], color='mediumspringgreen')
    
    # Plotting any potential test points
    if len(testPoints) != 0:
        ax.scatter(
        testPoints[:,0],testPoints[:,1],testPoints[:,2], color='black'
        )
    
    ax.legend()
    plt.show()

    return


# %% Defining Hubble volume and cosmic string wake wedge 유

# Defining ambient grid and Hubble lattice (volume)
hubbleScale = 6
hubbleLattice = hubbleScale * np.array([
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
     ])

# Defining the wake wedge 
deficitAngle = np.pi/2
wakeLength = 1
wakeTipPoint = np.array([0.5,0.5,0.5])
wakeEndPoints = np.array([
    [wakeTipPoint[0]+np.cos(deficitAngle/2),wakeTipPoint[1]+np.sin(deficitAngle/2),wakeTipPoint[2]],
    [wakeTipPoint[0]+np.cos(deficitAngle/2),wakeTipPoint[1]-np.sin(deficitAngle/2),wakeTipPoint[2]]
    ])
projectedWakePoints = np.array([wakeTipPoint+[0,0,wakeLength],wakeEndPoints[0]+[0,0,wakeLength],wakeEndPoints[1]+[0,0,wakeLength]])
middlePointWake = np.array([wakeTipPoint[0]+(1/2)*np.cos(deficitAngle/2),wakeTipPoint[1],wakeTipPoint[2]+wakeLength/2])
wakeWedge = np.array([
    wakeTipPoint, wakeEndPoints[0], wakeEndPoints[1], projectedWakePoints[0], projectedWakePoints[1], projectedWakePoints[2]
])

# Plotting Hubble lattice and wake wedge
quickPlot(
    hubbleScale,[hubbleLattice, wakeWedge],['Hubble lattice','Wake wedge'],
    ['blue','red']
    )

# Plotting again with smaller axis
quickPlot(
    hubbleScale/2,[hubbleLattice, wakeWedge, middlePointWake],['Hubble lattice','Wake wedge','Wake center'],
    ['blue','red','black']
    )


# %% Generating rotations via group representations and translations on wake wedge state 유

# Real representation of SO(3) group
def xAxisRotation(theta):
    return np.array([
        [1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]
    ])

def yAxisRotation(theta):
    return np.array([
        [np.cos(theta),0,np.sin(theta)],[0,1,0],[-np.sin(theta),0,np.cos(theta)]
    ])

def zAxisRotation(theta):
    return np.array([
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
    hubbleScale/2,[rotatedWedge, wakeWedge, middlePointWake],['Rotated wake','Wake wedge','Wake center'],
    ['darkorchid','red','black']
    )


#%% Computing the boundary of the wedge via convex hulls 유

# Using cube as an initial example (unscaled Hubble lattice)
cubeVertices = np.array([
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
     ])

# Plotting cube complex closure
quickPlotHull(cubeVertices)

# Repeating the above calculation with the wake wedge
quickPlotHull(wakeWedge)

    
#%% Checking if a set of points lie within the hull 유

def pointChecker(hullPoints, testPoints, plot=False):
    
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
        quickPlotHull(hullPoints,testPoints)
    
    else:
        
        return criterion
    

# Running preliminary examples
pointChecker(wakeWedge,np.array([[1,1,1]]),True)
pointChecker(wakeWedge,np.array([[2,1,1]]),True)
pointChecker(wakeWedge,np.array([[1,1,1],[2,1,1]]))


#%% Computing brightness temperature function δT(z) for points within wake 유

# Defining temporary function for brightness temperature
def brightnessTemperature(z):
    return np.sqrt(1+z)

# Finding grid of points within (unscaled) Hubble lattice (3 x 1000 elements)
meshGrid = np.mgrid[0:2:0.2, 0:2:0.2, 0:2:0.2]
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

for i in range(len(physicalPositions)):
    if pointChecker(wakeWedge,np.array([physicalPositions[i]])) == True:
        temperatureValues[i] = brightnessTemperature(physicalPositions[i][2])
        wedgePoints.append(physicalPositions[i])
        wedgeTemperatures.append(brightnessTemperature(physicalPositions[i][2]))
        
physicalWedgePoints = np.array(wedgePoints)
physicalWedgeTemperatures = np.array(wedgeTemperatures)

# Plotting physical points within wake excluding ambient points
quickPlot(
    hubbleScale/4,physicalWedgePoints,['Internal temperature points'],
    ['red'], physicalWedgeTemperatures
    )


# %% CODE GRAVEYARD: Testing rotations: A = dot(A, R.T) where A is an array of points, R.T is the transposed rotation matrix R
#region
xUnitVector = np.array([1,0,0])
yUnitVector = np.array([0,1,0])
zUnitVector = np.array([0,0,1])
rotationMatrix = xAxisRotation(-np.pi/2) #CW rotation about x-axis
rotatedXVector = np.dot(xUnitVector,rotationMatrix.T)
rotatedYVector = np.dot(yUnitVector,rotationMatrix.T) # Note: Rotations are CCW!
rotatedZVector = np.dot(zUnitVector,rotationMatrix.T)
print('Rotated x-direction unit vector:' , rotatedXVector) #gives expected (1,0,0)
print('Rotated y-direction unit vector:' , rotatedYVector) #gives expected (0,0,-1)
print('Rotated z-direction unit vector:' , rotatedZVector) #gives expected (0,1,0)

#Finding middle point of wake
middlePointWake = np.array([wakeTipPoint[0]+(1/2)*np.cos(deficitAngle/2),wakeTipPoint[1],wakeTipPoint[2]+wakeLength/2])

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_ylim(0,hubbleScale/2)
ax.set_xlim(0,hubbleScale/2)
ax.set_zlim(0,hubbleScale/2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.scatter(
    hubbleLattice[:,0],hubbleLattice[:,1],hubbleLattice[:,2], color='blue',
    label = 'Hubble lattice'
    )
ax.scatter(
    middlePointWake[0],middlePointWake[1],middlePointWake[2], color='green',
    label = 'Wake center'
    )
ax.scatter(
    wakeWedge[:,0],wakeWedge[:,1],wakeWedge[:,2],
    label = 'Wake wedge', color='red'
    )
ax.legend()


def modifiedPointChecker(points, p, plot=True):
    
    # Defining hull and test points
    hull = ConvexHull(points)
    new_points = np.append(points, p, axis=0)
    new_hull = ConvexHull(new_points)
    
    # Checking if points are within the hull
    if list(hull.vertices) == list(new_hull.vertices):
        criterion = True
    else:
        criterion = False
    
    if plot == True:
    
        # Plotting the hull & test points
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.plot(points.T[0], points.T[1], points.T[2], "ko", color='blue')
        
        for s in hull.simplices:
            s = np.append(s, s[0])
            ax.plot(points[s, 0], points[s, 1], points[s, 2], color='mediumspringgreen')
            
        for i in ["x", "y", "z"]:
            eval("ax.set_{:s}label('{:s}')".format(i, i))
            
        ax.scatter(
        p[:,0],p[:,1],p[:,2], color='black' #POTENTIALLY ADD PLOTTING BOOLEAN
        )
    
        print(criterion)
        return criterion
        
    else:
    
        print(criterion)
        
        return criterion

modifiedPointChecker(wakeWedge,np.array([[1,1,1]]),True)
modifiedPointChecker(wakeWedge,np.array([[2,1,1]]),False)
modifiedPointChecker(wakeWedge,np.array([[1,1,1],[2,1,1]]),False)

#Testing plotting & color maps
import matplotlib.pyplot as plt
import numpy as np

x = np.random.rand(5)
y = np.random.rand(5)

plt.scatter(x, y, c=[0,1,2,3,4])
plt.show()

#endregion