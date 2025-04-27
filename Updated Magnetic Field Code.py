#Jacob Lucian Underwood
#Magnetic Field Loop Animation Coding Project
#Rev: 12/26/2024

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

    #Pre-allocated Magnetic Field Components
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

    B = B.reshape(X.shape + (3,))
    Bx = B[...,0]
    By = B[...,1]
    Bz = B[...,2]
    MagField = [Bx, By, Bz]

    AnimateMagneticField(Radius,Segments,Thickness,MagField,save_animation=False)

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
    radSec =np.column_stack((Rad*np.cos(theta), -Rad*np.sin(theta),np.zeros_like(theta)))

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

def AnimateMagneticField(Rad,precision,thick,B,save_animation):

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

    ##Step 3: Creating the Shape of the Magnetic Field (Quiver Plot)

    #Had assistance from GPT to learn how to make the mesh and the quiver.

    #Creating the Magnetic Field Locations
    FieldX, FieldY, FieldZ = np.meshgrid(np.linspace(-1.5*Rad,1.5*Rad,10),
                                         np.linspace(-1.5*Rad,1.5*Rad,10),
                                         np.linspace(-1.5*Rad,1.5*Rad,10))

    #Adjusting Vector for easier readability
    Bx = B[0]
    By = B[1]
    Bz = B[2]

    #Plotting the Field
    quiver = ax.quiver(FieldX, FieldY, FieldZ, Bx, By, Bz,
                            alpha = .5, length = Rad/5, normalize = True)

    #Step 4: Animative the quiver. Help from GPT.
    def update(frame):
        nonlocal quiver
        quiver.remove()
        quiver = ax.quiver(FieldX, FieldY, FieldZ, Bx*np.sin(frame),By,Bz,
                               alpha = 0.5, length = Rad/4)
        #ax.view_init(elev = 30, azim = frame)
        return quiver,

    ani = FuncAnimation(fig,update,frames = np.arange(0,360,2), interval = 50)
    if save_animation:
        ani.save("magnetic_field.mp4", writer = "ffmpeg")
    
    
    plt.show()

main()
