#Jacob Lucian Underwood
#Magnetic Field Loop Animation Coding Project
#Rev: 12/25/2024



#Necessary Packages
import math #Original package. Probably will be removed due to utilizing numpy.
import numpy as np #A major component in making efficient code. Used for matrices/arrays.
from numpy.linalg import norm, cross #Needed for calculation.
import matplotlib.pyplot as plt #A major plotting component. Needed for models.
from matplotlib.animation import FuncAnimation #Needed for animating the code.

#The layout of the code:
def main():

    #Initial Conditions and Variables from Problem:
    Current = 4         #Current flowing through the Wire. (Amps)
    Radius = 0.04       #Radius of the Wire-Loop. (Meters)
    Thickness = 0.005   #The Thickness of the Wire. (Meters)
    Segments = 50       #How many parts the Wire is divided into for calculations.


    ##Points for Calculations:
    #Pre-Computing the Mesh:
    x_Range = np.linspace(-2*Radius, 2*Radius, 10)
    y_Range = np.linspace(-2*Radius, 2*Radius, 10)
    z_Range = np.linspace(-2*Radius, 2*Radius, 10)
    Points, X, Y, Z = ObservationPoints(x_Range,y_Range,z_Range)

    #Pre-allocationd Magnetic Field Components
    B = np.zeros((len(Points),3))

    #Biot-Savart Constants
    #These are placed here to avoid re-computation with each loop
    dSComp = -2*math.pi*Radius/Segments     #Constants for dS Calculations
    mu0 = 4*math.pi*10**(-7)                #Tesla*Meters/Seconds
    CrossComp = mu0*Current/(4*math.pi)     #dS-rHat Cross Product constants

    #The formatting for this section was partially-suggested by GPT. Originally it had
    #a different loop that make use of an enumerate command; as a whole I didn't
    #understand what it had written so I sought to deviate and at least still needed
    #assistance to make the code usable.
    for i in range(len(Points)):

        B[i,:] = (MagFieldCalc(Current,Points[i],Radius,Segments,dSComp,mu0,CrossComp))

    MagFieldPlot(Radius,Segments,Thickness,B)
    #MagField = B

    

def ObservationPoints(x,y,z):
    #Rather than using a singular point and determining the Magnetic Field there, I've
    #chosen to expand to the a set surrounding about the Loop. Hence the number of ranges.

    
    X,Y,Z = np.meshgrid(x,y,z)
    
    #Because of how python works, it works out easier to use the flatten command and
    #convert the large meshgrids into a single array.
    points = np.stack((X.flatten(),Y.flatten(),Z.flatten()),axis = -1)

    return points, X, Y, Z

def MagFieldCalc(I, Point, Rad,segments,dSComp,mu0,CrossComp):
    
    ###Initializing Various Arrays/Vectors and Constants
    #Angle between horizon and partition location
    theta = np.linspace(0,2*math.pi,segments)

    #Segments are effectively in 2D
    radSec = np.column_stack((Rad*np.cos(theta), -Rad*np.sin(theta),np.zeros_like(theta)))

    #Line Differentials wrt Segment Locations
    dS = np.column_stack((Rad*np.cos(theta),-Rad*np.sin(theta),np.zeros_like(theta)))

    #Distance Vector between Point of Interest and Segment Location
    rVector = Point-radSec

    #Magnitude of the Distance Vector
    rMagnitude = norm(rVector)

    #Distance Unit Vector
    rHat = rVector/rMagnitude

    #Magnetic Field Differential
    dB = cross(dS,rHat)*CrossComp/rMagnitude**2

    #Summing the components to receive Magnetic Field
    B = np.sum(dB,axis = 0)*1e6

    return B

def MagFieldPlot(Rad,precision,thick,B):

    ##Step 1: Creating the Shape of the Wire.

    #The Wire is in a loop; it's shape is that of a Torus (aka a glorified donut).

    #Code copied from mattiagiuri.com due to unfamiliarity of the parametric equations
    #for the shape
    
    #U and V are respective circular angles.
    #U is the anglular position of the Wire's segment.
    U = np.linspace(0,2*math.pi,precision)
    #V is the angle of the segment's circle.
    V = np.linspace(0,2*math.pi,precision)
    U,V = np.meshgrid(U,V)

    #Parametric Equations of the Wire
    WireX = (Rad + thick*np.cos(V))*np.cos(U)
    WireY = (Rad + thick*np.cos(V))*np.sin(U)
    WireZ = thick*np.sin(V)

    #Plotting the Wire
    fig = plt.figure()
    #Original Code had a mistake for the add_subplot command. Changed wrt to
    #Matplotlib tutorials.
    ax = fig.add_subplot(projection='3d')
    ax.plot_surface(WireX, WireY, WireZ, antialiased=True,color = 'gray',alpha = 0.5)
    
    ##Step 2: Creating the "Shape" of the Current (Just a Circle).
    #Parametric Equations of the Current
    CurrentX = Rad*np.cos(U)
    CurrentY = Rad*np.sin(U)
    CurrentZ = np.zeros_like((CurrentX)) #Current Loop is effectively 2-D

    #Setting up the gradient field
    Current_dX = np.gradient(CurrentX)
    Current_dY = np.gradient(CurrentY)
    Current_dZ = np.gradient(CurrentZ)
    #Plotting the Current
    ax.quiver(CurrentX, CurrentY, CurrentZ,
              Current_dX, Current_dY, Current_dZ,
              color = 'red')

    ##Step 3: Creating the Shape of the Magnetic Field (Coil wrapping around Torus)
    #Assistance from ChatGPT for coding the Coil

    #Choosing Precision to both be the Number of Coils and the Points per Coil.
    #Like U and V from Step 1, U is the angle for the whole Torus, V is for the
    #the segment's angle.
    theta = np.linspace(0,2*math.pi,precision**2)
    phi = np.linspace(0,2*math.pi*precision,precision**2)

    #Adjusting the Width of the Coil to be wider than the Wire
    thick  += thick/2

    #Parametric Equations of the Coil
    #CoilX = (Rad+thick*np.cos(phi))*np.cos(theta)
    #CoilY = (Rad+thick*np.cos(phi))*np.sin(theta)
    #CoilZ = thick*np.sin(phi)

    #Coil_dX = np.gradient(CoilX)
    #Coil_dY = np.gradient(CoilY)
    #Coil_dZ = np.gradient(CoilZ)
    #ax.quiver(CoilX,CoilY,CoilZ,Coil_dX,Coil_dY,Coil_dZ,color = 'blue',length = Rad*8)

    #Plotting the Coil
    #ax.plot(CoilX, CoilY, CoilZ, color = 'blue')

    #Creating the Magnetic Field
    FieldX, FieldY, FieldZ = np.meshgrid(np.linspace(-2*Rad,2*Rad,10),
                                         np.linspace(-2*Rad,2*Rad,10),
                                         np.linspace(-2*Rad,2*Rad,10))

    Bx = B[0,:]
    By = B[1,:]
    Bz = B[2,:]

    #Plotting the Field
    ax.quiver(FieldX, FieldY, FieldZ, Bx, By, Bz,
              alpha = .5, length = Rad/4, normalize = True)
    
    plt.show()

main()
