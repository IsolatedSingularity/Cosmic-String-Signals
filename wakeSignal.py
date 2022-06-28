#%% Importing modules & notes 유
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
# Make a plotting function + close sections of code
# Make fixed values parameters
# Testing convex hull 
# Connect lines & define boundary + new center of wake to old center
# Boolean function (boundary-volume criterion)
# Wake length dimensions
# Scale dimensions from Mpc to redshift (non-natural, astrophysical units) (take amplitude in 21cmfast and convert to natural untis)
# Combine: Resolution, frequency (21cmFAST Mpc units, output amplitude in Kelvin) & match grid points (EXPLICIT), Gmu free param (start with 10^-8)
# Add variational brightness temperature 
# Analytically how does crossover redshift depend on Gmu
# Add anisotropies on string manifold
# Clean up all code and plots
# Optimizations of space

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
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_ylim(0,hubbleScale)
ax.set_xlim(0,hubbleScale)
ax.set_zlim(0,hubbleScale)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.scatter(
    hubbleLattice[:,0],hubbleLattice[:,1],hubbleLattice[:,2], color='blue',
    label = 'Hubble lattice'
    )
ax.scatter(
    wakeWedge[:,0],wakeWedge[:,1],wakeWedge[:,2],
    label = 'Wake wedge', color='red'
    )
ax.legend()

# Plotting with smaller axis
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
    wakeWedge[:,0],wakeWedge[:,1],wakeWedge[:,2],
    label = 'Wake wedge', color='red'
    )
ax.scatter(
    middlePointWake[0],middlePointWake[1],middlePointWake[2],
    label = 'Wake center', color='black'
    )
ax.legend()


# %% Generating rotation via group representations and translations on wake wedge state 유

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
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.set_ylim(0,hubbleScale/2)
ax.set_xlim(0,hubbleScale/2)
ax.set_zlim(0,hubbleScale/2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.scatter(
    middlePointWake[0],middlePointWake[1],middlePointWake[2],
    label = 'Wake center', color='black'
    )
ax.scatter(
    rotatedWedge[:,0],rotatedWedge[:,1],rotatedWedge[:,2], color='darkorchid',
    label = 'Rotated wake'
    )
ax.scatter(
    wakeWedge[:,0],wakeWedge[:,1],wakeWedge[:,2], color='red',
    label = 'Unrotated wake', cmap='hot'
    )
# ax.scatter(
#     shiftedWedge[:,0],shiftedWedge[:,1],shiftedWedge[:,2], color='blue',
#     label = 'Shifted wake', cmap='hot'
#     )
ax.legend()


#%% Computing the boundary of the wedge via convex hulls 유

# Using cube as an initial example (unscaled Hubble lattice)
cubeVertices = np.array([
    [0,0,0],[1,0,0],[1,1,0],[0,1,0],
    [0,0,1],[1,0,1],[1,1,1],[0,1,1]
     ])

# Defining convex closure
hull = ConvexHull(cubeVertices)

# Creating subplots
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

# Plotting the cube vertices
ax.plot(cubeVertices.T[0], cubeVertices.T[1], cubeVertices.T[2], "ko", color='blue')

# Plotting simplices connecting to vertices (2 simplices per square face)
for s in hull.simplices:
    s = np.append(s, s[0])  # Permuting back to first coordinate
    ax.plot(cubeVertices[s, 0], cubeVertices[s, 1], cubeVertices[s, 2], color='mediumspringgreen')

# Defining plot axes labels
for i in ["x", "y", "z"]:
    eval("ax.set_{:s}label('{:s}')".format(i, i))

plt.show()


# Repeating the above calculation with the wake wedge
pts = wakeWedge
hull = ConvexHull(pts)

fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

ax.plot(pts.T[0], pts.T[1], pts.T[2], "ko", color='blue')

for s in hull.simplices:
    s = np.append(s, s[0])
    ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], color='mediumspringgreen')

for i in ["x", "y", "z"]:
    eval("ax.set_{:s}label('{:s}')".format(i, i))

plt.show()

    
#%% Checking if a set of points lie within the hull 유

def pointChecker(points, p):
    
    # Defining hull and test points
    hull = ConvexHull(points)
    new_points = np.append(points, p, axis=0)
    new_hull = ConvexHull(new_points)
    
    # Checking if points are within the hull
    if list(hull.vertices) == list(new_hull.vertices):
        criterion = True
    else:
        criterion = False
    
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
    p[:,0],p[:,1],p[:,2], color='black'
    )
    
    print(criterion)
    
    return criterion

# Running preliminary examples
pointChecker(wakeWedge,np.array([[1,1,1]]))
pointChecker(wakeWedge,np.array([[2,1,1]]))
pointChecker(wakeWedge,np.array([[1,1,1],[2,1,1]]))


#%% Computing brightness temperature function δT(z) for points within wake

# Defining temporary function for brightness temperature
def brightnessTemperature(z):
    return np.sqrt(1+z)

# Finding grid of points within (unscaled) Hubble lattice (3 x 1000 elements)
X, Y, Z = np.mgrid[0:1:0.1, 0:1:0.1, 0:1:0.1]
positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])

# Working with smaller grid for testing (3 x 8 elements)
X, Y, Z = np.mgrid[0:1:0.5, 0:1:0.5, 0:1:0.5]
positions = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])

# Getting array of all points within the space
emptyPositions = []
for i in range(len(positions.T)): #8
    emptyPositions.append([positions[0][i],positions[1][i],positions[2][i]])

physicalPositions = np.array(emptyPositions)

# Assigning brightness temperature to points within the wedge
temperatureValues = np.zeros((8,1))
for i in range(len(physicalPositions)):
    if pointChecker(wakeWedge,physicalPositions[i]) == True:
        temperatureValues[i] = brightnessTemperature(physicalPositions[i][2]) #Dimension ERROR









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
#endregion