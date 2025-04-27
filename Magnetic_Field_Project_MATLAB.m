%Section 1: Setting up the Variables

Current = (input('Enter the current going through your loop in Amps: '))*100 %Current running through the loop in amps
segments = 50 %the number of segments of loop needed for Riemann Sum
PointLocation = [0 0 0]
for i = 1:length(PointLocation)
    PointLocation(i) = (input('Enter origin (center of the loop) point:'))*100;
end
radius = (input('Enter the radius of your loop, in meters: '))*100 %radius of the loop in meters
mu0 = ((4*pi)*(10^(-7))) %Tesla meters per second

%Section 2: Getting the Components Ready
i = (0:segments/49:segments);

theta = 2*pi/segments*i; %Angle made as segments more
radComps = [radius*cos(theta); -radius*sin(theta); 0*i];
lineDifferential = [-2*pi*radius/segments*sin(theta); -2*pi*radius/segments*cos(theta); 0*i];


%Section 3: Calculating the r-Vector and r-Magnitude

rVec = PointLocation'-radComps;

for j =1:length(rVec)
    rMagnitude(j) = norm(rVec(:,j));
end

%Section 4: Calculating Magnetic Field differentials

%Below are the respective components to the Mag-Field Differentials
%Given that each individual piece is going to have differentials, we need to make the
%r-Hat components and Cross-Product Pieces into lists.

rHat = rVec./rMagnitude;
Crossed = cross(lineDifferential,rHat);
dB = mu0*Current*Crossed./(4*pi*(rMagnitude.^2));

Magnetic_Field_Result = sum(dB')*10^6

%Section 5: Animating the current flowing through  the loop

%First, setting up the shape.
%The shape will be of a perfectly circular loop of the original radius.
%Its crossections will be a simple point.

%Finding the values for the Wire Loop
xWireLoopComp = radius*cos(theta);
yWireLoopComp = radius*sin(theta);

%Finding the values for the magnetic field loops
xMagLoopComp = xWireLoopComp+1;
yMagLoopComp = yWireLoopComp+1;
zMagLoopComp = ones([1 50]);

%Plotting the wire loop
figure

plot3(xWireLoopComp,yWireLoopComp,0*(1:segments))
hold on

%Plotting the segments at which the magnetic field will be shown
for i = 1:50
    plot3(xWireLoopComp(i),yWireLoopComp(i),0*i,'r.',MarkerSize=10)
    hold on
    for l = 1:50
        plot3(xMagLoopComp(l),yMagLoopComp(l),zMagLoopComp(l),'g-',MarkerSize = 10)
    end
    drawnow
end