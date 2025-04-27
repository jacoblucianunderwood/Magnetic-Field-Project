#Jacob Lucian Underwood
#Magnetic Field Loop Animation Coding Project
#Rev: 4/24/2025, Fixing up the code for simplicity.

##Program Purpose:
#Originally, this program was created as a bonus project for Murray State's PHY 255 course.
#The original inception date was the Spring for 2022 as that's when I took PHY 255.
#The goal at the time was to calculate the magnetic field at a determined point away from
#a current carrying loop. This was done by using Biot-Savart's Law and summing up the total
#contribution from the entire loop to garner the Magnetic Field at that point.
#Not long afterwards, I had a decided to improve upon this program in many ways, but notably
#through adding animations and in general making the code more efficient (The orignal code
#did not make good use of functions).

##What has changed:
#The first notable change is the amount of packages utilized. Before, I was only allowed to
#use the math package. Now we also tack on numpy and mpl.
#The second notable change is the use of functions. Back then I didn't want to use them.
#I was just pretty lazy and thought I would clean up the code later. Later ended up being
#nearly 3 years later.
#Lastly, I've added code to create a model showing the loop, the current, and the magnetic
#field.

##Goals:
#For now, the goal is to simplify the code further and also include animations. Once the
#animations have been incorportated that should be sufficient.

#Necessary Packages
import numpy as np #A major component in making efficient code. Used for matrices/arrays.
from numpy.linalg import norm, cross #Needed for calculation.
import matplotlib.pyplot as plt #A major plotting component. Needed for models.

#The layout of the code:
def main():

    #Initial Conditions and Variables from Problem:
    Current = 4         #Current flowing through the Wire. (Amps)
    Radius = 0.04       #Radius of the Wire-Loop. (Meters)
    Thickness = 0.005   #The Thickness of the Wire. (Meters)
    Segments = 50       #How many parts the Wire is divided into for calculations.


    ##Points for Calculations:
    #Pre-Computing the Mesh:
    Range = np.linspace(-2*Radius, 2*Radius, 4) #Computing for a cube area, so all 3 lengths are the same
    Points = ObservationPoints(Range)
    #Points = np.array([[0.04,0.02,0.01]])


    #Pre-allocated Magnetic Field Components
    B = np.zeros(Points.shape)

    #Biot-Savart Constants
    #These are placed here to avoid re-computation with each loop
    dSComp = -2*np.pi*Radius/Segments     #Constants for dS Calculations
    mu0 = 4*np.pi*10**(-7)                #Tesla*Meters/Seconds
    CrossComp = mu0*Current/(4*np.pi)     #dS-rHat Cross Product constants

    #Attempt to Vectorize Later?
    B = MagFieldCalc(Current,Points,Radius,Segments,dSComp,mu0,CrossComp)

    #for i in range(len(Points)):
        #print(Points[i])

     #   B[i,:] = (MagFieldCalc(Current,Points[i],Radius,Segments,dSComp,mu0,CrossComp))
        #print(B)
        #B[0,:]: X values
        #B[1,:]: Y Values
        #B[2,:]: Z Values

    MagneticFieldPlot(Radius,Segments,Thickness,B[:,0],B[:,1],B[:,2])


def ObservationPoints(Range):
    #Rather than using a singular point and determining the Magnetic Field there, I've
    #chosen to expand to the a set surrounding about the Loop. Hence the number of ranges.
    X, Y, Z = np.meshgrid(Range,Range,Range)
    
    #Because of how python works, it works out easier to use the flatten command and
    #convert the large meshgrids into a single array.
    #Had assistance from GPT for this
    points = np.vstack((X.ravel(),Y.ravel(),Z.ravel())).T

    return points

def MagFieldCalc(I, Point, Rad,segments,dSComp,mu0,CrossComp):
    
    ###Initializing Various Arrays/Vectors and Constants
    #Angle between horizon and partition location
    theta = np.linspace(0,2*np.pi,segments,endpoint = False)

    #Segments are effectively in 2D
    radSec = np.column_stack((Rad*np.cos(theta), -Rad*np.sin(theta),np.zeros_like(theta)))

    #Line Differentials wrt Segment Locations
    dS = np.column_stack((-2*Rad*np.pi/segments*np.sin(theta),-2*Rad*np.pi/segments*np.cos(theta),np.zeros_like(theta))) 
    
    if Point.shape[0]>1:

        #Reshaping the matrices for rVector:
        PointDims = Point.shape
        Point = Point.reshape(1, PointDims[1],PointDims[0])
        radSec = np.tile(radSec[:,:,np.newaxis],(1,1,PointDims[0]))
        #Reshaping the matrix for the later crossproduct:
        dS = np.tile(dS[:,:,np.newaxis],(1,1,PointDims[0]))

    

    #Distance Vector between Point of Interest and Segment Location
    rVector = Point-radSec

    #Magnitude of the Distance Vector
    rMagnitude = norm(rVector,axis = 1,keepdims=True)

    #Distance Unit Vector
    rHat = rVector/rMagnitude

    

    #Magnetic Field Differential
    dB = cross(dS,rHat, axis = 1)*CrossComp/rMagnitude**2

    #Summing the components to receive Magnetic Field
    B = np.sum(dB,axis = 0)*1e6
    #print(B)

    return B

def MagneticFieldPlot(Rad,precision,thick,Bx,By,Bz):

    ##Step 1: Creating the Shape of the Wire.

    #The Wire is in a loop; it's shape is that of a Torus (aka a glorified donut).

    #Code copied from mattiagiuri.com due to unfamiliarity of the parametric equations
    #for the shape
    
    #U and V are respective circular angles.
    #U is the anglular position of the Wire's segment.
    U = np.linspace(0,2*np.pi,precision)
    #V is the angle of the segment's circle.
    V = np.linspace(0,2*np.pi,precision)
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

    ##Step 3: Creating the Shape of the Magnetic Field (Quiver Plot)

    #Had assistance from GPT to learn how to make the mesh and the quiver.

    #Creating the Magnetic Field Locations
    Field = np.linspace(-2*Rad,2*Rad,len(Bx)) #Creating a shape like area
    FieldX,FieldY,FieldZ = np.meshgrid(Field,Field,Field)
    print(len(Bx))

    #Plotting the Field
    step = 1
    quiver = ax.quiver( FieldX[::step, ::step, ::step,],
                        FieldY[::step, ::step, ::step],
                        FieldZ[::step, ::step, ::step],
                        Bx[::step], By[::step], Bz[::step],
                        alpha = .5, length = Rad/5, normalize = True)
    
    
    plt.show()

main()
