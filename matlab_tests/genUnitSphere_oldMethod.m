%{

The formula for a point on the sphere is as follows:
P = ( 
    sin(phi) cos(theta), 
    sin(phi) sin(theta),
    cos(theta) 
    )

%}

%computes the interval values, method 1
%interval = pi/30;
%angleValues = 0:interval:2*pi;
%angleValuesHalf = 0:interval:pi/2;


numValues = 50;

%computes random theta and phi values, method 2
angleValues = rand(1,numValues)*(2*pi);
angleValuesHalf = rand(1,numValues)*(pi/2);

%makes the meshgrid for method 1 and 2
%[Phi, Theta] = meshgrid(angleValuesHalf,angleValues);

%makes random phi,theta pair, method 3
Phi = rand(1,numValues*numValues)*(2*pi);
Theta = rand(1,numValues*numValues)*(2*pi);

%computes the sphere
pointsX = sin(Phi).*cos(Theta);
pointsY = sin(Phi).*sin(Theta);
pointsZ = cos(Phi);

plot3(pointsX,pointsY,pointsZ,'r.');