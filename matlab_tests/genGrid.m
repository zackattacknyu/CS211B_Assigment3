numVals = 10;
interval = 2/numVals;

xVals = -1:interval:1;
yVals = -1:interval:1;

[X,Y] = meshgrid(xVals,yVals);

%draw a circle
Theta = 0:pi/100:2*pi;
Xcircle = cos(Theta);
Ycircle = sin(Theta);

hold on
plot(X,Y,'b');
plot(Y,X,'b');
plot(Xcircle,Ycircle,'g','LineWidth',3);

for iteration = 1:50
    randomValsX = rand(size(X)).*interval + X;
    randomValsY = rand(size(Y)).*interval + Y;
    
    squaredDist = randomValsX.^2 + randomValsY.^2;
    
    randomValsX = randomValsX(squaredDist <= 1);
    randomValsY = randomValsY(squaredDist <= 1);

    plot(randomValsX,randomValsY,'r.');
end

hold off