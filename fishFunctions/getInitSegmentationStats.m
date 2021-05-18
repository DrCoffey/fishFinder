%Computes region precision (overlap) and boundary precision 
function [REC, PR, F1, Overlap,BP,BDE,RECL,PRL,F1L,LGT] = getInitSegmentationStats(GT,Seg,LSeg)
lines = size(GT,1);
cols = size(GT,2);
GT0 = GT;
GT = convertGTBoundarytoRegions( GT );
CC = bwconncomp(GT, 4);
LGT = labelmatrix(CC);
SS = regionprops(LGT, 'Area');
GT = ismember(LGT, find([SS.Area] >= 100));
CC = bwconncomp(GT, 4);
LGT = labelmatrix(CC);

c = 0;
w = 0;
for i=1:lines
    for j=1:cols,
        if GT(i,j) == 1 && Seg(i,j) == 1,
            c = c+1;
        elseif GT(i,j) == 1 || Seg(i,j) == 1,
            w = w+1;
        end
    end
end

REC = c/sum(GT(:));
PR = c/sum(Seg(:));
F1 = 2*REC*PR/(REC+PR);
Overlap = c / (c+w); %overlap
            
[BDE, ~,BP] = compare_image_boundary_error(double(GT), double(Seg));

if ~isempty(LSeg),
    SGT = regionprops(LGT, 'Centroid');
    SS = regionprops(LSeg, 'Centroid');
    Detlabels = [];
    for i=1:length(SGT),
        u = round(SGT(i).Centroid(1));
        v = round(SGT(i).Centroid(2));
        if LSeg(v,u) > 0,
           Detlabels = [Detlabels LSeg(v,u)];
        else
          %  disp(sprintf('(%d,%d)-%d',v,u,LGT(v,u)));
        end
    end
    Detlabels = unique(Detlabels);
    
    TP = length(Detlabels);
    FP = length(SS)-TP;
    FN = length(SGT)-TP;
    RECL = TP/(TP+FN);
    PRL = TP/(TP+FP);
    F1L = 2*RECL*PRL/(RECL+PRL);
end

end



%A MATLAB Toolbox
%
%Compare two segmentation results using the Boundary Displacement Error
%
%IMPORTANT: The input two images must have the same size!
%
%Authors: John Wright, and Allen Y. Yang
%Contact: Allen Y. Yang <yang@eecs.berkeley.edu>
%
%(c) Copyright. University of California, Berkeley. 2007.
%
%Notice: The packages should NOT be used for any commercial purposes
%without direct consent of their author(s). The authors are not responsible
%for any potential property loss or damage caused directly or indirectly by the usage of the software.

function [averageError, returnStatus,BP] = compare_image_boundary_error(imageLabels1, imageLabels2)

returnStatus = 0;
[imageX, imageY] = size(imageLabels1);
if imageX~=size(imageLabels2,1) || imageY~=size(imageLabels2,2)
    error('The sizes of the two comparing images must be the same.');
end

if isempty(find(imageLabels1~=imageLabels1(1), 1))
    % imageLabels1 only has one group
    boundary1 = zeros(size(imageLabels1));
    boundary1(1,:) = 1;
    boundary1(:,1) = 1;
    boundary1(end,:) = 1;
    boundary1(:,end) = 1;
else
    % Generate boundary maps
    [cx,cy] = gradient(imageLabels1);
    [boundaryPixelX{1},boundaryPixelY{1}] = find((abs(cx)+abs(cy))~=0);
    
    boundary1 = abs(cx) + abs(cy) > 0;
end

if isempty(find(imageLabels2~=imageLabels2(1), 1))
    % imageLabels2 only has one group
    boundary2 = zeros(size(imageLabels2));
    boundary2(1,:) = 1;
    boundary2(:,1) = 1;
    boundary2(end,:) = 1;
    boundary2(:,end) = 1;    
else    
    % Generate boundary maps
    [cx,cy] = gradient(imageLabels2);
    [boundaryPixelX{2},boundaryPixelY{2}] = find((abs(cx)+abs(cy))~=0);
    
    boundary2 = abs(cx) + abs(cy) > 0;
end

% boundary1 and boundary2 are now binary boundary masks. compute their
% distance transforms:
D1 = bwdist(boundary1);
D2 = bwdist(boundary2);
% compute the distance of the pixels in boundary1 to the nearest pixel in
% boundary2:
dist_12 = sum(sum(boundary1 .* D2 ));
dist_21 = sum(sum(boundary2 .* D1 ));

avgError_12 = dist_12 / sum(sum(boundary1));
avgError_21 = dist_21 / sum(sum(boundary2));

averageError = (avgError_12 + avgError_21) / 2;
BP = (sum(sum(boundary1)) + sum(sum(boundary2)))/(dist_12+dist_21);
end

