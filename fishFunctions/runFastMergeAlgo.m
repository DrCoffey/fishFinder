
function [EL,BestIClust,BestTotalPerf,NUMEllipses,area,AIC,BIC,minAICBIC,SI,nCompl,hfig1] = runFastMergeAlgo(A,lines,cols,AICBIC_SELECTION,Overlaping)
%disp('FAST MERGING');
[nCompl] = getObjectComplexity(A);
INF = 10000000000;
Cent = [];
Rad = [];
hfig1 = [];
%SK = bwmorph(A,'skel',inf);
SK = bwmorph(A,'thin',inf);%????
EP = bwmorph(SK,'endpoints');
toPlot = 0;
BD = SK.*bwdist(1-A);
FULLPLOT = 0;
stats = regionprops(double(A), 'Area');
area = stats.Area;

[Points(:,1) Points(:,2)] = find(BD > 0);
W = zeros(1,size(Points,1));
for i=1:length(W),
    W(i) = BD(Points(i,1),Points(i,2));
end

%LOCAL MAXIMA COMPUTATION
LM = zeros(1,size(Points,1));
VAL = zeros(1,size(Points,1));
ite = 0;
for i=1:size(Points,1),
    v = [(Points(:,1)-Points(i,1)).^2+(Points(:,2)-Points(i,2)).^2];
    maxDist = sqrt(max(v)/(mean(v)+0.001));
    
    ite = ite+1;
    LM(ite) = i;
    VAL(ite) = BD(Points(i,1),Points(i,2))+(1/(1+1000*maxDist));
end
LM = LM(1:ite);
VAL = VAL(1:ite);

[SVAL,pos] = sort(VAL,'descend');
ite = 1;
i = 1;
Cent(ite,1:2) = Points(LM(pos(i)),1:2);
Rad(ite) = SVAL(i);
%Selection of candicate circle centers
for i=2:length(VAL),
    v = sqrt([(Cent(:,1)-Points(LM(pos(i)),1)).^2+(Cent(:,2)-Points(LM(pos(i)),2)).^2]);
    if SVAL(i) < 0.03*SVAL(1)+1,
        break;
    end
    ok = 1;
    x0 = Points(LM(pos(i)),1);
    y0 = Points(LM(pos(i)),2);
    
    for j=1:length(v),
        mdist = v(j);
        if Overlaping < 0.5
            if EP(x0,y0) == 0 && mdist <= Rad(j)+SVAL(i)%+SVAL(i),%min(Rad(j),SVAL(i))+1%Rad(j)+SVAL(i)+1, %out of cirles; %max(Rad(j),SVAL(i))+1,  %
                ok = 0;
                break;
            elseif EP(x0,y0) == 1 && mdist <= Rad(j),
                ok = 0;
                break;
            end
        else
            if EP(x0,y0) == 0 && mdist <= Rad(j)%+SVAL(i),%min(Rad(j),SVAL(i))+1%Rad(j)+SVAL(i)+1, %out of cirles; %max(Rad(j),SVAL(i))+1,  %
                ok = 0;
                break;
            elseif EP(x0,y0) == 1 && mdist <= Rad(j)/2,
                ok = 0;
                break;
            end
        end
    end
    if ok == 1,
        ite = ite+1;
        Cent(ite,1:2) = Points(LM(pos(i)),1:2);
        Rad(ite) = SVAL(i);
    end
    
end

[EL0,~] = initEll(Cent,Rad,[1:length(Rad)]);
[EL,~,~,~] = runEllClustering(EL0,A,area);
[EL,IClust,DTemp,TotalPerf] = runEllClustering(EL,A,area);

BestDTemp = DTemp;
BestEL = EL;
ELLSET = [1:length(Rad)];
BEST_ELLSET = ELLSET;
BETTERSOL = 0;
BestIClust = IClust;
BestTotalPerf = TotalPerf;
[AIC(length(ELLSET)),BIC(length(ELLSET)),RES_AICBIC,bestAICBIC,SI(length(ELLSET),1:3)] = getAIC_BIC(nCompl,TotalPerf,length(ELLSET),AICBIC_SELECTION,IClust,EL,ELLSET);
minAICBIC = RES_AICBIC;
if getTotalPerfMultfromArea(EL,1:length(ELLSET),A) == 1
    firstSol = 1;
else
    firstSol = 0;
end

if toPlot == 1,
    [hfig1] = drawDistEllClusteting(BD(lines+1:2*lines,cols+1:2*cols)+double(A(lines+1:2*lines,cols+1:2*cols)),EL0,lines,cols);
    title(sprintf(' Complexity = %2.2f  NUMBER OF ELLIPSES = %d',nCompl,length(ELLSET)));
end

if FULLPLOT == 1,
    [~] = drawDistEllClusteting(BestDTemp(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
    title(sprintf('TotalPerf = %2.4f - %d',TotalPerf,length(BEST_ELLSET)));
end

for ite=1:length(Rad)-1,  
    [ConnMatrix] = getConnMatrix(IClust,max(ELLSET));
    M = size(ConnMatrix,1);
    AICBICgain = INF*ones(M,M);
    newEL = [];
    
    for i=1:M,
        for j=i+1:M,
            lab1 = i;
            lab2 = j;
            if EL(i).Label == i && EL(j).Label == j% && ConnMatrix(lab1,lab2) > 0, %check for merging
                [COVERgain,newEL] = checkMergeEll(A,newEL,IClust,lab1,lab2,EL,DTemp);
                TotalPerfAfter = (COVERgain+TotalPerf*area)/ area;
                AICBICgain(lab1,lab2) = getAICBICgain(nCompl,TotalPerf,TotalPerfAfter,AICBIC_SELECTION);
            end
        end
    end
%     
    change = 0;

    %merging process
   
    temp = min(min(AICBICgain));
    if temp >= INF,
        continue;
    end
    merges = 0;
    mergeSet = [];
    while 1,
        temp = min(min(AICBICgain));
        
        if temp <= 0 || merges == 0,
            [lab1,lab2] = find(AICBICgain == temp);
            lab1 = lab1(1);
            lab2 = lab2(1);
            
            merges = merges+1;
            AICBICgain(lab1,:) = INF;
            AICBICgain(:,lab1) = INF;
            AICBICgain(lab2,:) = INF;
            AICBICgain(:,lab2) = INF;
            ELLSET = setdiff(ELLSET, lab2);
            EL(lab2).Label = lab1;
            EL(lab2).ELLSET = [];
            for i=lab2:numel(EL),
                if ~isempty(EL(i).ELLSET),
                    EL(i).ELLSET = EL(i).ELLSET-1;
                end
            end
            mergeSet = [mergeSet lab1 lab2];
            EL(lab1) = newEL(lab1,lab2);
            IClust(IClust == lab2) = lab1;
        else
            break;
        end
        
    end
    mergeSet;
    NUMEllipses = length(ELLSET);
    
    [EL,IClust,DTemp,TotalPerf] = runEllClusteringForMerge(EL,ELLSET,IClust,area);
    %TotalPerf = sum([EL.InArea]) / area;
    [AIC1,BIC1,RES_AICBIC,bestAICBIC,SI1(1,1:3)] = getAIC_BIC(nCompl,TotalPerf,NUMEllipses,AICBIC_SELECTION,IClust,EL,ELLSET);
    
    AIC(length(ELLSET)) = AIC1;
    BIC(length(ELLSET)) = BIC1;
    SI(length(ELLSET),1:3) = SI1(1,1:3);
    
    [gtp,overlap] = getTotalPerfMultfromArea(EL,ELLSET,A);
    
    if (RES_AICBIC < minAICBIC && gtp == 1) || firstSol == 0,
        if gtp == 1,
            firstSol = 1;
        end
        BestIClust = IClust;
        BEST_ELLSET = ELLSET;
        minAICBIC = RES_AICBIC;
        change = 1;
        BETTERSOL = 0;
        BestDTemp = DTemp;
        BestTotalPerf = TotalPerf;
        BestEL = EL;
        if FULLPLOT == 1,
            [~] = drawDistEllClustetingMerge(BestDTemp(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
            title(sprintf('TotalPerf = %2.4f - %d',BestTotalPerf,length(BEST_ELLSET)));
        end
    else
        BETTERSOL = BETTERSOL+1;
        BestDTemp = DTemp;
    end
    
    
    if numel(EL) == 1
%         BETTERSOL
%         bestAICBIC
%         minAICBIC
        break;
    end
    
    if change == 0,%?????10
        if FULLPLOT == 1,
            [~] = drawDistEllClustetingMerge(BestDTemp(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
            title(sprintf('(Retry %d), TotalPerf = %2.4f - %d',BETTERSOL,TotalPerf,length(ELLSET)));
        end
    end
end
[EL] = informEL(BestEL);
NUMEllipses = numel(EL);
end


%Initilization of Ellipses
function [COVERgain,newEL] = checkMergeEll(BW0,newEL,IClust,lab1,lab2,EL,BestDTemp)
lines = size(IClust,1);
cols = size(IClust,2);
BW = (IClust == lab1 | IClust == lab2);

if sum(sum(BW)) == 0,
    COVERgain = -size(IClust,1)*size(IClust,2);
    return;
end
% 
% [x y] = meshgrid(1:max(lines,cols),1:max(lines,cols));
% X0 = EL(lab1).C(1);
% Y0 = EL(lab1).C(2);
% el=((x-X0)/EL(lab1).a).^2+((y-Y0)/EL(lab1).b).^2<=1;
% el = rotateAround(el,Y0,X0,EL(lab1).phi,'nearest');
% el1 = el(1:lines,1:cols);
% X0 = EL(lab2).C(1);
% Y0 = EL(lab2).C(2);
% el=((x-X0)/EL(lab2).a).^2+((y-Y0)/EL(lab2).b).^2<=1;
% el = rotateAround(el,Y0,X0,EL(lab2).phi,'nearest');
% el2 = el(1:lines,1:cols);
% el = max(el1,el2);
% el = min(el,BW0);
% BW = max(el,BW);
stats = regionprops(double(BW), 'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation');


area = stats.Area;
C = stats.Centroid;
e =  stats.MajorAxisLength / stats.MinorAxisLength;
X0 = C(1);
Y0 = C(2);
phi = stats.Orientation;


a = sqrt(e*area/pi);
b = a/e;

[x y] = meshgrid(1:max(lines,cols),1:max(lines,cols));

el=((x-X0)/a).^2+((y-Y0)/b).^2<=1;

el = rotateAround(el,Y0,X0,phi,'nearest');
el = el(1:lines,1:cols);

p1 = [];
p2 = [];
[p1(:,1) p1(:,2)] = find(el == 1 & BW == 1);
[p2(:,1) ~] = find(el == 1 | BW == 1);

tomh_area = size(p1,1) / area;
tomh_enwsh = size(p1,1) / size(p2,1);

BWELL = BW == 1 & BestDTemp <= 1;
Coverage_bef = sum(sum(BWELL));
Coverage_after = size(p1,1);
COVERgain = Coverage_after-Coverage_bef;

newEL(lab1,lab2).a = a;
newEL(lab1,lab2).b = b;
newEL(lab1,lab2).C = C;
newEL(lab1,lab2).phi = phi;
newEL(lab1,lab2).InArea = size(p1,1);
newEL(lab1,lab2).outPixels = size(p2,1) - size(p1,1);
newEL(lab1,lab2).tomh_area = tomh_area;
newEL(lab1,lab2).tomh_enwsh = tomh_enwsh;
newEL(lab1,lab2).Label = lab1;
newEL(lab1,lab2).ELLSET = EL(lab1).ELLSET;

end


%Matrix of connected regions 
function [ConnMatrix] = getConnMatrix(IClust,M)

lines = size(IClust,1);
cols = size(IClust,2);

ConnMatrix = zeros(M,M);

[x,y] = find(IClust > 0);

for i=1:length(x),
    a = x(i);
    b = y(i);
    lab1 = IClust(a,b);
    for i1=-1:2:1,
        for i2=-1:2:1,
            u = a+i1;
            v = b+i2;
            if u >= 1 && v >= 1 && u <= lines && v <= cols,
                if IClust(u,v) > 0 && IClust(a,b) ~= IClust(u,v),
                    lab2 = IClust(u,v);
                    ConnMatrix(lab1,lab2) = ConnMatrix(lab1,lab2)+1;
                    ConnMatrix(lab2,lab1) = ConnMatrix(lab2,lab1)+1;
                end
            end
        end
    end
end
end


function [AICBICgain] = getAICBICgain(nCompl,TotalPerfBef,TotalPerfAfter,AICBIC_SELECTION)

CONST = 1;
MODEL_PAR = 1;

AICBICgain = nCompl*log(1-TotalPerfAfter)-nCompl*log(1-TotalPerfBef)-MODEL_PAR*CONST*log(nCompl);

if AICBIC_SELECTION == 1,
    AICBICgain = nCompl*log(1-TotalPerfAfter)-nCompl*log(1-TotalPerfBef)-2*CONST*MODEL_PAR;
end
end


%update the EL set after the merging of lab2
function [EL] = informEL(BestEL)
k = 1;
for i=1:numel(BestEL),
    if BestEL(i).Label == i,
        EL(k) = BestEL(i);
        k = k+1;
    end
end

end

