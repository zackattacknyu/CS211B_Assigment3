%{

The formula for a point on the sphere is as follows:
P = ( 
    sin(phi) cos(theta), 
    sin(phi) sin(theta),
    cos(theta) 
    )

%}
clear all;
numGridVals = 30;
numIterations = 5;

interval = 2/numGridVals;
intervalValues = -1:interval:1;

sizeIntVals = size(intervalValues);
numValues = sizeIntVals(2);

%{
The random number generator generates from a normal distribution
   where the standard deviation is standardDev and the mean is mean. 
This way the numbers tend toward the center of the box
%}
standardDev = 1/8;
mean = 1/2;

%makes the meshgrid
[gridXvals, gridYvals] = meshgrid(intervalValues,intervalValues);

%the normal used to get the hemisphere we want
normal = [1/sqrt(2),1/sqrt(2),0];
%normal = [0,0,1];

aVal = normal(1);
bVal = normal(2);
cVal = normal(3);

for iteration = 1:numIterations
    
    randomIntervalValues = rand(1,numValues).*interval;

    %do uniform distribution
    Xvals = gridXvals + rand(numValues,numValues).*interval;
    Yvals = gridYvals + rand(numValues,numValues).*interval;

    %do normal distribution
    %Xvals = gridXvals + (randn(numValues,numValues).*standardDev + mean).*interval;
    %Yvals = gridYvals + (randn(numValues,numValues).*standardDev + mean).*interval;

    %only put it in the center
    %Xvals = gridXvals + 0.5.*interval;
    %Yvals = gridYvals + 0.5.*interval;

    radius = 1;
    squaredDist = Xvals.^2 + Yvals.^2;
    XvalsPlot = Xvals(squaredDist <= radius^2);
    YvalsPlot = Yvals(squaredDist <= radius^2);
    ZvalsPlot = sqrt(1 - XvalsPlot.^2 - YvalsPlot.^2);
    
    XvalsPlot = [XvalsPlot XvalsPlot];
    YvalsPlot = [YvalsPlot YvalsPlot];
    ZvalsPlot = [ZvalsPlot -ZvalsPlot];
    
    %{
    To rotate the points, all points on the unit sphere are specified, 
        and to do only the ones with respect to the normal, I used only
        the vectors (a,b,c) where (a,b,c)(x,y,z) >= 0. 
    %}
    dotProds = XvalsPlot.*aVal + YvalsPlot.*bVal + ZvalsPlot.*cVal;
    
    XvalsToPlot = XvalsPlot(dotProds >= 0);
    YvalsToPlot = YvalsPlot(dotProds >= 0);
    ZvalsToPlot = ZvalsPlot(dotProds >= 0);
    

    hold on
    plot3(XvalsToPlot,YvalsToPlot,ZvalsToPlot,'r.');
    hold off
    
end

