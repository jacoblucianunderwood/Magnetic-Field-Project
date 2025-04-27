#Jacob Lucian Underwood
#Magnetic Field Loop Animation Coding Project
#12/22/2024

import math

#Section 1: Setting up the Variables

I = 4 
segments = 50 #the number of segments of loop needed for Riemann Sum
point = [0.04, 0.02, 0.01]
radius = 0.04 
mu0 = ((4*math.pi)*(10**(-7))) #Tesla meters per second

#Section 2: Getting the Components Ready

radSec = [] #List of sections of circle. Each item of the lsit is the coordinates
dS = [] #List of line differential components (List of a list)
theta = [] #Angle made as segments more

for i in range(segments):
    theta.append(2*math.pi/segments*i)
    #theta is in radians. given there are a number of segments, we split the whole angle
    #into 50 pieces and multiply by i for which individual segment you're on
    radSec.append([radius*math.cos(theta[i]),-radius*math.sin(theta[i]),0])
    dSComp = -2*math.pi*radius/segments
    dS.append([dSComp*math.sin(theta[i]), dSComp*math.cos(theta[i]),0])
    # for purposes of code, line differential always in X-Y

#print(theta)
#Section 3: Calculating the r-Vector and r-Magnitude

rVec = [] #List of X-Components for the r-Vector
rMag = [] #List of Magnitudes of r-Vector
for i in range(segments):
    rVec.append([point[0]-radSec[i][0],point[1]-radSec[i][1],point[2]])
    #Change in X-Components to find the r-Vector's X-Component
    rMag.append(math.sqrt((rVec[i][0])**2+(rVec[i][1])**2+(rVec[i][2])**2))

#Section 4: Calculating Magnetic Field differentials

#Below are the respective components to the Mag-Field Differentials
#Given that each individual piece is going to have differentials, we need to make the
#r-Hat components and Cross-Product Pieces into lists.

dB = [] #List of Magnetic Field Differentials
rHat = [] #rHat Values
MagConst = []
Cross = []

for i in range(segments):
    #First, we need to set up the rHat components in order to compute the cross product pieces
    rHat.append([rVec[i][0]/rMag[i],(rVec[i][1])/rMag[i],rVec[i][2]/rMag[i]])
    #Now, we compute the cross product pieces
    Cross.append([dS[i][1]*rHat[i][2]-dS[i][2]*rHat[i][1],
                  dS[i][2]*rHat[i][0]-dS[i][0]*rHat[i][2],
                  dS[i][0]*rHat[i][1]-dS[i][1]*rHat[i][0]])
    
    #Finally, we compute the Magnetic Field Differential Components
    MagConst.append(mu0*I/(4*math.pi*rMag[i]**2))
    dB.append([Cross[i][0]*MagConst[i], Cross[i][1]*MagConst[i],Cross[i][2]*MagConst[i]])

print(dS[20])
    
#Section 5: Summing up the Magnetic Field Differential Components to find the Magnetic Field Vector

B = [0]*3

for i in range(segments):
    B[0] += dB[i][0]*10**6
    B[1] += dB[i][1]*10**6
    B[2] += dB[i][2]*10**6

#print('The magnetic Field should be: ', B)
