%Computes region precision (overlap) and boundary precision
function [Jaccard, MAD, Hausdorff, DiceFP,DiceFN,FP,FN,LGT] = getStats(GT,Seg,LSeg)
lines = size(GT,1);
cols = size(GT,2);
GT0 = GT;
se = strel('disk', 1);
%Labels
GT = convertGTBoundarytoRegions( GT );
GTSEG = max(GT,GT0);
GTSEG = imerode(GTSEG,se);

CC = bwconncomp(GT, 4);
LGT = labelmatrix(CC);

SS = regionprops(LGT, 'Area');
GT = ismember(LGT, find([SS.Area] >= 120)); %ground truth for splitting
CC = bwconncomp(GT, 4);
LGT = labelmatrix(CC); %ground truth for splitting (LABELS)

CC = bwconncomp(GTSEG, 4);
GT = labelmatrix(CC);
SS = regionprops(GT, 'Area');
GTSEG = ismember(GT, find([SS.Area] >= 120)); %ground truth for segmentation

Max1 = max(GTSEG,Seg);
Min1 = min(GTSEG,Seg);
tomh = sum(Min1(:));
enosh = sum(Max1(:));

Jaccard = 100*tomh / enosh;
R_ = 1-GTSEG;
ar = min(R_,Seg);
ar = sum(ar(:));
par = GTSEG+Seg;
par = sum(par(:));
DiceFP = 2*100*ar/par;

R_ = 1-Seg;
ar = min(GTSEG,R_);
ar = sum(ar(:));
DiceFN = 2*100*ar/par;
hgt = zeros(1,max(LGT(:)));
hgtE = hgt;
hseg = zeros(1,max(LSeg(:)));
hsegE = hseg;
for i=1:lines,
    for j=1:cols,
        l = LGT(i,j);
        if l > 0,
            hgt(l) = hgt(l)+Min1(i,j);
            hgtE(l) = hgtE(l)+Max1(i,j);
        end
        l = LSeg(i,j);
        if l > 0,
            hseg(l) = hseg(l)+Min1(i,j);
            hsegE(l) = hsegE(l)+Max1(i,j);
        end
    end
end
hgt = hgt./max(1,hgtE);
hseg = hseg./max(1,hsegE);

tempSeg = Seg;

for i=1:lines,
    for j=1:cols,
        l = LGT(i,j);
        if l > 0 && hgt(l) < 0.5,
           % tempGT(i,j) = 0;
            tempSeg(i,j) = 0;
        end
        l = LSeg(i,j);
        if l > 0 && hseg(l) < 0.5,
          %  tempGT(i,j) = 0;
            tempSeg(i,j) = 0;
        end
    end
end

[Gmag] = edge(GTSEG,'Sobel',0.001);
Gmag = Gmag > 0;
%Gmag = imerode(Gmag,se);
[P1, P2] = find(Gmag >0);
P = zeros(length(P1),2);
P(:,1) = P1;
P(:,2) = P2;

[Gmag] = edge(tempSeg,'Sobel',0.001);
Gmag = Gmag > 0;
%Gmag = imerode(Gmag,se);
[P1, P2] = find(Gmag >0);
Q = zeros(length(P1),2);
Q(:,1) = P1;
Q(:,2) = P2;
CC = bwconncomp(imdilate(Seg,strel('disk', 2)), 4);
LISEG = labelmatrix(CC);

[Hausdorff,MAD] = HausdorffMADDist(Q,P,LISEG);%????


FN = 0;
FP = 0;
%SPLITTING - EVALUATION
if ~isempty(LSeg),
    CC = bwconncomp(Seg, 4);
    LISEG = labelmatrix(CC);

    CC = bwconncomp(GTSEG, 4);
    LIGT = labelmatrix(CC);
    
    SSeg = regionprops(LSeg, 'Centroid');
    SGT = regionprops(LGT, 'Centroid');
    SGTSEG = regionprops(LISEG, 'Centroid');
    
    pseg = zeros(length(SSeg),5);
    for i=1:length(SSeg),
        if isnan(SSeg(i).Centroid(1)) == 0,
            pseg(i,1:2) = round(SSeg(i).Centroid);
            pseg(i,3) = LISEG(pseg(i,2),pseg(i,1));
            pseg(i,4) = LIGT(pseg(i,2),pseg(i,1));
            pseg(i,5) = LGT(pseg(i,2),pseg(i,1));
        end
    end
    
    pgt = zeros(length(SGT),5);
    for i=1:length(SGT),
        if isnan(SGT(i).Centroid(1)) == 0,
            pgt(i,1:2) = round(SGT(i).Centroid);
            pgt(i,3) = LISEG(pgt(i,2),pgt(i,1));
        end
    end
    
    numDetectedObjects = zeros(1,length(SGTSEG));
    numGTObjects = zeros(1,length(SGTSEG));
    for i=1:length(SSeg),
        if pseg(i,3) > 0 && pseg(i,4) > 0,
            id = pseg(i,3);
            numDetectedObjects(id) = numDetectedObjects(id)+1;
        end
    end
    for i=1:length(SGT),
        if pgt(i,3) > 0,
            id = pgt(i,3);
            numGTObjects(id) = numGTObjects(id)+1;
        end
    end
    FP = find(numGTObjects >= 1 & numDetectedObjects >= 1 & numDetectedObjects > numGTObjects);
    FP = length(FP);
    FN = find(numGTObjects >= 1 & numDetectedObjects >= 1 & numDetectedObjects < numGTObjects);
    FN = length(FN);
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

