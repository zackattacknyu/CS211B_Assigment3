%{

The formula for a point on the sphere is as follows:
P = ( 
    sin(phi) cos(theta), 
    sin(phi) sin(theta),
    cos(theta) 
    )

%}

%computes the interval values
%interval = pi/30;
%angleValues = 0:interval:2*pi;
%angleValuesHalf = 0:interval:pi/2;

numValues = 30;
angleValues = rand(1,30)*(2*pi);
angleValuesHalf = rand(1,30)*(pi/2);

%makes the meshgrid
[Phi, Theta] = meshgrid(angleValuesHalf,angleValues);

%computes the sphere
pointsX = sin(Phi).*cos(Theta);
pointsY = sin(Phi).*sin(Theta);
pointsZ = cos(Phi);

plot3(pointsX,pointsY,pointsZ,'r.');