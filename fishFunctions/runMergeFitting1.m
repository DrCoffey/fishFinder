%It runs the DEFA (Merging)method
%I: original binary image
%AICBIC_SELECTION: Set AICBIC_SELECTION = 1, to use AIC is selected else BIC is used
%IClust:Clustering of image pixels
%EL: Ellipses 
%NUMEllipses:Number of Ellipses 

function [IClust,EL,NUMEllipses] = runMergeFitting1(I,AICBIC_SELECTION)
Iorig = I;
lines = size(I,1);
cols = size(I,2);
I = zeros(3*lines,3*cols);
I(lines+1:2*lines,cols+1:2*cols) = Iorig;
toPlot = 0;
IClust = I;

[EL,IClust,TotalPerf,NUMEllipses,area] = runFastMergeAlgo1(IClust);

if toPlot == 1,
    [ok] = drawEllClusteting(IClust(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
    title(sprintf('%4.2f',100*TotalPerf));


    drawEllClusteting(IClust(lines+1:2*lines,cols+1:2*cols),EL,lines,cols);
    for i=1:NUMEllipses,
        hold on;
        text(EL(i).C(1) -cols,EL(i).C(2)-lines,sprintf('%d',i));
    end
    title(sprintf('%4.2f',100*TotalPerf));

    figure;
    plot(AIC,'-.o');
    hold on;
    plot(BIC,'--rs');
    legend('AIC','BIC');
end

IClust = IClust(lines+1:2*lines,cols+1:2*cols);


end

function [EL,BestIClust,BestTotalPerf,NUMEllipses,area] = runFastMergeAlgo1(A)
%disp('FAST MERGING');
[nCompl] = getObjectComplexity(A);
Cent = [];
Rad = [];
%SK = bwmorph(A,'skel',inf);
SK = bwmorph(A,'thin',inf);%????
BD = SK.*bwdist(1-A);
stats = regionprops(double(A), 'Area');
area = stats.Area;
BestIClust = A;
BestTotalPerf = 1;
NUMEllipses = 1;

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


[EL0,~] = initEll(Cent,Rad,[1]);
[EL,~,~,~] = runEllClustering(EL0,A,area);
[EL,IClust,DTemp,TotalPerf] = runEllClustering(EL,A,area);

end


