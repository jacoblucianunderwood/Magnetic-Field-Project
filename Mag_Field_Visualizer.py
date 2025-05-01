#Jacob Lucian Underwood
#Magnetic Field Visualizer Project
#Rev: 5/01/2025: Rewriting the Code to find issues.

#Code Purpose:
#Originally this code was designed as a simple calculator for a Physics 2 Bonus Project on Biot-Savart's Law

#Now it has become an expanded project focused on visualizing the Magnetic Field about a current-carrying wire loop.

#Current Challenges:
#1. Getting a good visualization of the Magnetic Field. But given the fact every element of the circuit is producing 
#   a loop, getting a visualization may be difficult
#2. Vectorizing the code and getting accurate calculations performed.

#Necessary Libraries and functions for simplification
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import cross, norm

def main():
    
    ##Initializing Conditions
    Current = 4         #Amps (Current through Wire)
    Radius = 0.04       #Meters (Radius of Wire Loop)
    Thickness = 0.005   #Meters (Thickness of Wire)
    Segments = 50       #Number of Wire Segments for Calculations


    ##Points for Calculations
    #Pre-Computing the Observation Volume (Cubic Volume):
    Range = np.linspace(-2*Radius,2*Radius,15)#Cube Dimensions
    Points=ObservationPoints(Range)

    #Test Case to Verify Math
    #Points = np.array([0.04,0.02,0.01])

    ##Magnetic Field Computations
    #Constants:
    dSComp = -2*np.pi*Radius/Segments      #Constants for dS Calculations
    mu0 = 4*np.pi*10**(-7)                 #Tesla-Meters/Second (Magnetic Constant)
    CrossComponent = mu0*Current/(4*np.pi) #dS-rHat Cross Product Constanta

    #Computation:
    B_Vectors = MagFieldCalc(Current,Points,Radius,Segments,dSComp,mu0,CrossComponent)

    ##Plotting:
    MagFieldPlot(B_Vectors,Points,Segments,Thickness,Radius)
    



def ObservationPoints(Range):

    #This function is used to gather a wide range of 3-D points for computations

    #Creating the Cubic Volume
    X,Y,Z = np.meshgrid(Range,Range,Range)

    #"Flattening" the dimensions for calculations later
    points = np.vstack((X.flatten(),Y.flatten(),Z.flatten())).T

    return points

def MagFieldCalc(I,Point,Radius,Segments,dSComp,mu0,CrossComp):

    ##Initializing Vectors and Constants

    #Angular positions of Segments
    theta = np.linspace(0,2*np.pi,Segments,endpoint = False) 

    #Radius Coordinates of Segments
    rad_Segment = Radius*np.column_stack((np.cos(theta),-np.sin(theta),np.zeros_like(theta)))
    
    #Line Differentials of Segments
    dS = -2*Radius*np.pi/Segments*np.column_stack((-np.sin(theta),-np.cos(theta),np.zeros_like(theta)))

    if Point.shape[0] > 1:
        #Reshaping the matrices for rVector and Cross Products:
        PointDims = Point.shape
        Point = Point.reshape(1,PointDims[1],PointDims[0])
        rad_Segment = np.tile(rad_Segment[:,:,np.newaxis],(1,1,PointDims[0]))
        dS = np.tile(dS[:,:,np.newaxis],(1,1,PointDims[0]))
 
    #Distance Vector between Point of Interest and Segment Location
    rVector = Point-rad_Segment

    #Magnitude of rVector
    #We choose Axis = 1 here to specifically target the 3D Components
    rMagnitude = norm(rVector,axis = 1, keepdims = True)

    #rVector Unit Vector
    rHat = rVector/rMagnitude

    #Magnetic Field Differential
    
    #We choose Axis = 1 here to specifically target the 3D Components
    dB = cross(dS,rHat,axis = 1)*CrossComp/rMagnitude**2

    #Total Magnetic Field at Points:
    B = np.sum(dB,axis = 0)*1e6

    return B

def MagFieldPlot(B_Vec,Point,Segment,Thick,Rad):

    ##First, we plot the Wire:

    #Wire is in the shape of a donut. Both relevant angles are that
    #of a circle, so we can use u_1 for simplicity.
    u_1 = np.linspace(0,2*np.pi,Segment)
    U_1, V_1 = np.meshgrid(u_1,u_1)
    #U_1 is of the Wire's Angular position at Segment
    #V_1 is the angle of that particular segment's circle

    #Wire Parametric Equations:
    WireX = (Rad + Thick*np.cos(V_1))*np.cos(U_1)
    WireY = (Rad + Thick*np.cos(V_1))*np.sin(U_1)
    WireZ = Thick*np.sin(V_1)

    #Plotting the Wire
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    ax.plot_surface(WireX,WireY,WireZ,color = 'k')

    ##Next, we plot the Magnetic Field Lines as a Quiver
    #Quick Adjustment to B_Vec for dimensional homogeneity
    B_Vec = B_Vec.T

    skip = np.arange(0,Point.shape[0],2)

    ax.quiver(Point[skip,0],Point[skip,1],Point[skip,2],
    B_Vec[skip,0],B_Vec[skip,1],B_Vec[skip,2],color = 'b',alpha = 0.25,normalize=True,arrow_length_ratio = 0.3,length = 0.005)

    plt.show()






main()