%{

The formula for a point on the sphere is as follows:
P = ( 
    sin(phi) cos(theta), 
    sin(phi) sin(theta),
    cos(theta) 
    )

%}

clear all;
maxNumValues = 100;

errorAmountArray = zeros(maxNumValues-30);
numValsArray = zeros(maxNumValues-30);
index = 1;
for numGridVals = 30:maxNumValues
    
    interval = 2/numGridVals;
    intervalValues = -1:interval:1;

    sizeIntVals = size(intervalValues);
    numValues = sizeIntVals(2);

    %makes the meshgrid
    [gridXvals, gridYvals] = meshgrid(intervalValues,intervalValues);

    randomIntervalValues = rand(1,numValues).*interval;

    %do uniform distribution
    Xvals = gridXvals + 0.5*interval;
    Yvals = gridYvals + 0.5*interval;

    radius = 1;
    squaredDist = Xvals.^2 + Yvals.^2;
    XvalsPlot = Xvals(squaredDist <= radius^2);
    YvalsPlot = Yvals(squaredDist <= radius^2);
    ZvalsPlot = sqrt(1 - XvalsPlot.^2 - YvalsPlot.^2);
    %plot3(XvalsPlot,YvalsPlot,ZvalsPlot,'r.');

    functionValues = 1./ZvalsPlot;
    approxPatchAreas = functionValues.*(interval*interval);
    totalSurfaceArea = sum(approxPatchAreas);
    
    errorAmountArray(index) = abs(totalSurfaceArea-2*pi);
    numValsArray(index) = numGridVals;
    index = index + 1;
end

plot(numValsArray,errorAmountArray);


