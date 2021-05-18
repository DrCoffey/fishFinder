function [EL,area,TotalPerf] = getBestFitEllipses(I,EL,NUMEllipses,area)
p = [];
s = 0;
for k=1:NUMEllipses,
    [EL,~,p1] = getBestFitEllipse(I,EL,k);
    p = union(p,p1);
   % s = s+length(p1);
end
%TotalPerf = s/area;
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
