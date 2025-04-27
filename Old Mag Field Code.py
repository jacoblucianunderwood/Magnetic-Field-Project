import math

#Section 1: Setting up the Variables

I = 4 # float(input('Enter the current going through your loop in Amps: ')) #Current running through the loop in amps
segments = 50 #the number of segments of loop needed for Riemann Sum
#xComponent = float(input('Enter the distance from the origin (center of the loop) along the x-value for your point, in meters: '))
#yComponent = float(input('Enter the distance from the origin (center of the loop) along the y-value for your point, in meters: '))
#zComponent = float(input('Enter the distance from the origin (center of the loop) along the z-value for your point, in meters: '))
point = [0.04, 0.02, 0.01] #[xComponent,yComponent,zComponent] #location of point in cartesian
radius = 0.04 #float(input('Enter the radius of your loop, in meters: ')) #radius of the loop in meters
mu0 = ((4*math.pi)*(10**(-7))) #Tesla meters per second
#Section is correct

#Section 2: Getting the Components Ready

radX = [] #List of radius x-components
radY = [] #List of radius y-components
dsX = [] #List of change in x direction (x-line differential)
dsY = [] #List of change in y direction (y-line differential)
dsZ = [0]*50 #No change in Z-direction
theta = [] #Angle made as segments more

for i in range(segments):
    theta.append(2*math.pi/segments*i)
    #theta is in radians. given there are a number of segments, we split the whole angle
    #into 50 pieces and multiply by i for which individual segment you're on
    radX.append(radius*math.cos(theta[i])) #determining x-component at each interval
    radY.append(-radius*math.sin(theta[i])) #determining y-component at each interval
    dsX.append(-2*math.pi*radius/segments*math.sin(theta[i])) #determining x-comp of the line differential
    dsY.append(-2*math.pi*radius/segments*math.cos(theta[i])) #determining y-comp of the line differential

#Section 3: Calculating the r-Vector and r-Magnitude

rVecX = [] #List of X-Components for the r-Vector
rVecY = [] #List of Y-Components for the r-Vector
rVecZ = [] #List of Z-Components for the r-Vector
rMagnitude = [] #List of Magnitudes of r-Vector
for i in range(segments):
    rVecX.append(point[0]-radX[i]) #Change in X-Components to find the r-Vector's X-Component
    rVecY.append(point[1]-radY[i]) #Change in y-Components to find the r-Vector's Y-Component
    rVecZ.append(point[2]) #Change in z-Components to find the r-Vector's Z-Component
    rMagnitude.append(math.sqrt(rVecX[i]**2+rVecY[i]**2+rVecZ[i]**2))

#Section 4: Calculating Magnetic Field differentials

#Below are the respective components to the Mag-Field Differentials
#Given that each individual piece is going to have differentials, we need to make the
#r-Hat components and Cross-Product Pieces into lists.

dBx = [] #X-Component of Magnetic Field Differential
dBy = [] #Y-Component of Magnetic Field Differential
dBz = [] #Z-Component of Magnetic Field Differential
xHat = [] #X-Component of rHat
yHat = [] #Y-Component of rHat
zHat = [] #Z-Component of rHat
iCross = [] #iHat-Vector Product from Cross Product of dS Vector cross rHat
jCross = [] #jHat-Vector Product from Cross Product of dS Vector cross rHat
kCross = [] #kHat-Vector Product from Cross Product of dS Vector cross rHat

for i in range(segments):
    #First, we need to set up the rHat components in order to compute the cross product pieces
    xHat.append((rVecX[i])/rMagnitude[i])
    yHat.append((rVecY[i])/rMagnitude[i])
    zHat.append((rVecZ[i])/rMagnitude[i])
    #Now, we compute the cross product pieces
    iCross.append(dsY[i]*zHat[i]-dsZ[i]*yHat[i])
    jCross.append(dsZ[i]*xHat[i]-dsX[i]*zHat[i])
    kCross.append(dsX[i]*yHat[i]-dsY[i]*xHat[i])
    #Finally, we compute the Magnetic Field Differential Components
    dBx.append((mu0*I*iCross[i])/(4*math.pi*rMagnitude[i]**2))
    dBy.append((mu0*I*jCross[i])/(4*math.pi*rMagnitude[i]**2))
    dBz.append((mu0*I*kCross[i])/(4*math.pi*rMagnitude[i]**2))
    
#Section 5: Summing up the Magnetic Field Differential Components to find the Magnetic Field Vector

B = [0]*3
convert = 1/(10**(-6))
for i in range(segments):
    B[0] += dBx[i]*convert
    B[1] += dBy[i]*convert
    B[2] += dBz[i]*convert

print('The magnetic Field should be: ', B)
