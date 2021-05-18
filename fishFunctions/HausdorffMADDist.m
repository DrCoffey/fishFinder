function [hd,MAD] = HausdorffMADDist(P,Q,LSeg)
% Calculates the Hausdorff Distance between P and Q
%
% hd = HausdorffDist(P,Q)
% [hd D] = HausdorffDist(P,Q)
% [hd D] = HausdorffDist(...,lmf)
% [hd D] = HausdorffDist(...,[],'visualize')
%
% Calculates the Hausdorff Distance, hd, between two sets of points, P and
% Q (which could be two trajectories). Sets P and Q must be matrices with
% an equal number of columns (dimensions), though not necessarily an equal
% number of rows (observations).
%
% The Directional Hausdorff Distance (dhd) is defined as:
% dhd(P,Q) = max p c P [ min q c Q [ ||p-q|| ] ].
% Intuitively dhd finds the point p from the set P that is farthest from
% any point in Q and measures the distance from p to its nearest neighbor
% in Q.
% 
% The Hausdorff Distance is defined as max{dhd(P,Q),dhd(Q,P)}
%
% D is the matrix of distances where D(n,m) is the distance of the nth
% point in P from the mth point in Q.
%
% lmf: If the size of P and Q are very large, the matrix of distances
% between them, D, will be too large to store in memory. Therefore, the
% function will check your available memory and not build the D matrix if
% it will exceed your available memory and instead use a faster version of
% the code. If this occurs, D will be returned as the empty matrix. You may
% force the code to forgo the D matrix even for small P and Q by calling the
% function with the optional 3rd lmf variable set to 1. You may also force
% the function to return the D matrix by setting lmf to 0. lmf set to []
% allows the code to automatically choose which mode is appropriate.
%
% Including the 'vis' or 'visualize' option plots the input data of
% dimension 1, 2 or 3 if the small dataset algorithm is used.
%
% Performance Note: Including the lmf input increases the speed of the
% algorithm by avoiding the overhead associated with checking memory
% availability. For the lmf=0 case, this may represent a sizeable
% percentage of the entire run-time.
%
%

% %%% ZCD Oct 2009 %%%
% Edits ZCD: Added the matrix of distances as an output. Fixed bug that
%   would cause an error if one of the sets was a single point. Removed
%   excess calls to "size" and "length". - May 2010
% Edits ZCD: Allowed for comparisons of N-dimensions. - June 2010
% Edits ZCD: Added large matrix mode to avoid insufficient memory errors
%   and a user input to control this mode. - April 2012
% Edits ZCD: Using bsxfun rather than repmat in large matrix mode to
%   increase performance speeds. [update recommended by Roel H on MFX] -
%   October 2012
% Edits ZCD: Added a plotting function for visualization - October 2012
%

sP = size(P); sQ = size(Q);

if ~(sP(2)==sQ(2))
    error('Inputs P and Q must have the same number of columns')
end

maxP = zeros(1,max(LSeg(:)));           % initialize our max value
MAD = maxP;
Par = maxP;
% loop through all points in P looking for maxes
for p = 1:sP(1)
    % calculate the minimum distance from points in P to Q
    l = LSeg(P(p,1),P(p,2));
    minP = min(sum( bsxfun(@minus,P(p,:),Q).^2, 2));
    if minP>maxP(l)
        % we've discovered a new largest minimum for P
        maxP(l) = minP;
    end
    MAD(l) = MAD(l)+minP;
    Par(l) = Par(l)+1;
end
MAD = MAD./max(Par,1);

% repeat for points in Q
% maxQ = 0;
% for q = 1:sQ(1)
%     minQ = min(sum( bsxfun(@minus,Q(q,:),P).^2, 2));
%     if minQ>maxQ
%         maxQ = minQ;
%     end
% end
maxP = sqrt(maxP);
MAD = sqrt(MAD);
hd = mean(maxP(Par > 0));
MAD = mean(MAD(Par > 0));
    
