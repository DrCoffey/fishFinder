%areaLim = 250
function [ rate,overlapMax ] = getTotalPerfMultfromArea(EL,SET,I)
global Constraints 
areaLim = Constraints(1);
minArea = EL(SET(1)).InArea;
maxArea = 0;
for i=1:length(SET),
    k = SET(i);
    minArea = min(minArea,EL(k).InArea);
    maxArea = max(maxArea,EL(k).InArea);
end

[~, overlapMax ] = getOverlapRatio(EL,SET,I);

rate = 1;

if minArea < areaLim,
    rate = 0.01;%(minArea/areaLim)^2;
end

g = minArea/maxArea;

if g < Constraints(2),
    rate = 0.01;%rate*(g/0.1)^2;
end

if overlapMax > Constraints(3),
    rate = 0.01;%rate*(g/0.1)^2;
end


