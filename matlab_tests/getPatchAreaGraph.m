Phi = 0:pi/100:pi/2;
Phi2 = Phi+pi/100;

deltaTheta = pi/50;

Area = (cos(Phi2)-cos(Phi)).*(-deltaTheta);

plot(Phi,Area);