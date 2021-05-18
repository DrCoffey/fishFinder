% Computes the object complexity n based on skeleton
% n is a positve number that is related to the object complexity indipedent
% from object size and orientation
function [n] = getObjectComplexity(A)
toPlot = 0;
B = zeros(size(A,1)+6,size(A,2)+6);
B(4:size(B,1)-3,4:size(B,2)-3) = A;
A = B;
SK = bwmorph(A,'thin',inf);
BW = bwdist(1-A);
BD = SK.*BW;
AInit = A;

AV = mean(BD(BD > 0));
rad = double(ceil(0.5*AV));
rad = min(rad,3);

if rad >=1,
    se = strel('disk',rad);
    A = imclose(A,se);
end
rad;

% 

%SK = bwmorph(A,'skel',inf);
SK = bwmorph(A,'thin',inf);
BW = bwdist(1-A);
BD = SK.*BW;
WEIGH = [];
[Points(:,1) Points(:,2)] = find(BD > 0);
S = bwboundaries(A);
Sarea = sqrt(length(find(A > 0)));
Degree = [];
rc = 0;

for i=1:length(S),
    if length(S{i}) > 8,
        rc = rc + 1;
    end
    Cirles = rc-1;
end
Cirles = max(Cirles,0);
W = zeros(1,size(Points,1));
for i=1:length(W),
    W(i) = BD(Points(i,1),Points(i,2));
end
W0 = W;

LM = zeros(1,size(Points,1));
UNSET = [];
ite = 0;

BP = bwmorph(SK,'branchpoints');
EP = bwmorph(SK,'endpoints');

[x,y] = find(BP == 1);

for i=1:size(x,1),
    v = [(Points(:,1)-x(i)).^2+(Points(:,2)-y(i)).^2];
    g1 = find(v == 0);%[3 x 3] block
    g = find(v <= max(2+0.00*Sarea,0.00*W0(i)*W0(i)));%[3 x 3] block

    if  length(g1) > 0, %branchpoints
        
        if ~isempty(find(UNSET == g1(1), 1)),
            continue;
        end
        
        ite = ite+1;
        LM(ite) = g1(1);
        WEIGH(ite) = 1;
        Degree(ite) = length(g)-1;
        UNSET = union(UNSET,g1); 
        continue;
    end
end

[x,y] = find(EP == 1);

for i=1:size(x,1),
    v = [(Points(:,1)-x(i)).^2+(Points(:,2)-y(i)).^2];
    g1 = find(v == 0);%[3 x 3] block
    g = find(v <= max(2+0.00*Sarea,0*W0(i)*W0(i)));%[3 x 3] block
    
    if  length(g1) > 0, %endpoints
        if ~isempty(find(UNSET == g1(1), 1)),
            continue;
        end
        ite = ite+1;
        LM(ite) = g1(1);
        WEIGH(ite) = 1;
        Degree(ite) = 1;
        UNSET = union(UNSET,g1); 
        continue;
    end
end


LM = LM(1:ite);

n = sum(WEIGH) + 2*Cirles;
[n2] = getShannonComplexity(Points,LM,W0,Degree);


if toPlot == 1,
    cmap = colormap(jet(256));
    cmap(1,:) = [1 1 1];
    figure;
    imagesc(BD+double(AInit));
    colormap(cmap);
    
    
    for i=1:length(LM),
        u = LM(i);
        hold on;
        if Degree(i) > 1,
            text(Points(u,2),Points(u,1),sprintf('%d',i));
            %text(Points(u,2),Points(u,1),sprintf('%d-(D%d)',i,Degree(i)));
        else
            text(Points(u,2),Points(u,1),sprintf('%d',i));
        end
    end
    
    title(sprintf(' S = %2.2f   Graph size = %d - Cirles = %d Points = %d',n2,n,Cirles,ite));
end

n = max(n2,1);



function [com] = getShannonComplexity(Points,LM,W0,Degree)

A = zeros(size(Points,1),size(Points,1));
BIGDIST = 10*size(Points,1)+10;
%A is the connection matrix 
for i=1:size(Points,1),
    for j=i+1:size(Points,1),
        d = (Points(i,1)-Points(j,1))^2 +(Points(i,2)-Points(j,2))^2;
        if d <= 2,
            A(i,j) = 0.75*d*(d-1)+1;%+BIGDIST*(length(find(LM == i))+length(find(LM == j)));
            A(j,i) = 0.75*d*(d-1)+1;%+BIGDIST*(length(find(LM == i))+length(find(LM == j)));
        end
    end
end


SA = sparse(A);

[dist] = graphallshortestpaths(SA);

Label = zeros(size(Points,1),2);

if length(LM) > 2,
    Ind = zeros(1,size(Points,1));
    for i=1:size(Points,1),
        if length(find(LM == i)) > 0,
            [~,pos] = find(LM == i);
            Label(i,1:2) = pos;
            Ind(i) = pos;
            continue;
        end
        
        vec = zeros(1,length(LM));
        for j=1:length(LM),
            vec(j) = dist(i,LM(j));
        end
        [~, pos1] = min(vec);
        vec(pos1) = BIGDIST^2;
        [~, path] = graphshortestpath(SA, i);
        
        for j=1:length(LM),
            if length(find(path{LM(j)} == LM(pos1))) > 0,
                vec(j) = BIGDIST*(vec(j)+1);
            end
        end
        [~, pos2] = min(vec);
        Label(i,1) = min(pos1,pos2);
        Label(i,2) = max(pos1,pos2);
    end
    
    
    tempDeg = Degree;
    ite = size(Label,1);
    ite0 = ite;
    CLM = zeros(length(LM),length(LM));
    ok = 0;
    while sum(tempDeg) > 0,
        sum(tempDeg);
        ok = min(ok,0);
        for i=1:length(LM),
            grade = min(tempDeg(tempDeg > 0));
            if tempDeg(i) == grade,
                [d, path] = graphshortestpath(SA, LM(i));
                [~,pos] = sort(d);
                
                for j=1:length(pos),
                    k = Ind(pos(j));
                    if k > 0 && tempDeg(k) > 0 && CLM(i,k) == 0 && i ~= k,
                        ite = ite+1;
                        SA(LM(i),LM(k)) = BIGDIST;
                        SA(LM(i),LM(k)) = BIGDIST;
                    
                        Label(ite,1) = min(i,k);
                        Label(ite,2) =  max(i,k);
                        W0(ite) = W0(LM(i));
                        ite = ite+1;
                        Label(ite,1) =  min(i,k);
                        Label(ite,2) =  max(i,k);
                        W0(ite) = W0(LM(k));
                        tempDeg(k) = tempDeg(k)-1;
                        tempDeg(i) = tempDeg(i)-1;
                        CLM(i,k) = 0;
                        CLM(k,i) = 0;
                        ok = ok+1;
                        break;
                    end
                end
            end
        end
        ok = ok-1;
        if ok < -5,
            %tempDeg
            break;
        end
    end
%     ite0
%     sum(Degree)
%     ite
end

L1 = unique(Label(:,1));
L2 = unique(Label(:,2));
n2 = zeros(size(L1,1),size(L2,1));
Rmax = max(W0);
Rmin = min(W0);
N = 16;



for i=1:size(L1,1)
    for j=1:size(L2,1)
        vec = find(Label(:,1) == L1(i) & Label(:,2) == L2(j));
        if length(vec) > 1,
            W1 = W0(vec);
            W = floor(N*(W1-Rmin) / (0.00000000001+Rmax-Rmin))+1;
            Pr = zeros(1,N);
            for k=1:length(W),
                Pr(W(k)) = Pr(W(k))+1;
            end
            Pr = Pr/sum(Pr);
            
            n2(i,j) = -sum(Pr.*log2(Pr+0.00000000001));
        end
    end
end
com = sum(sum(n2))+log2(size(Points,1)+0.0000001);


