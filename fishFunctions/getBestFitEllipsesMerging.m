function [EL,area,TotalPerf] = getBestFitEllipsesMerging(I,EL,NUMEllipses,area)
p = [];
for k=1:NUMEllipses,
    [EL,~,p1] = getBestFitEllipse(I,EL,k);
    p = union(p,p1);
end
TotalPerf = size(p,1)/area;

%AREA CHECK
% minArea = area;
% maxArea = 0;
% for k=1:NUMEllipses,
%     minArea = min(minArea,EL(k).InArea);
%     maxArea = max(maxArea,EL(k).InArea);
% end
% 
% if minArea < 250,
%     TotalPerf = 0.01;
% end
% 
% if minArea/maxArea < 0.1,
%     TotalPerf = 0.01;
% end
