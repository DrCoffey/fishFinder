%Performs clustering based on ellipse set EL 
function [EL,IClustNew,Dtemp,TotalPerf] = runEllClustering(EL,IClust,area)

IClustNew = IClust;
lines = size(IClust,1);
cols = size(IClust,2);
NUMEllipses = numel(EL);
Dtemp = zeros(lines,cols);
ite = 0;
Thresh_D = 3;
while 1,
    ite = ite+1;
    changes = 0;
%     k1 = 0;
%     p = zeros(1,2*area);

    for i=lines/3:2*lines/3,
        for j=cols/3:2*cols/3,
            if IClust(i,j) > 0,
                d = zeros(1,NUMEllipses);
                for k=1:NUMEllipses,
                    OAdist = norm([j i] - EL(k).C);%?????
                     ration = OAdist/max(EL(k).a,0.00001);
                     
                    if  ration > Thresh_D 
                        OXdist = (EL(k).a+EL(k).b)/2;
                    else
                        OXdist = getOX([j i],EL(k));
                    end
                    %OXdist = getOX([j i],EL(k));
                    d(k) = OAdist/max(OXdist,0.00001);
%                     if d(k) <= 1,
%                         k1 = k1+1;
%                         p(k1) = j+i*lines*cols;
%                     end
                end
                [Dtemp(i,j), pos] = min(d);
               
                if pos ~= IClustNew(i,j),
                   changes = changes+1; 
                end
                IClustNew(i,j) = pos;
            end
        end
    end
    Thresh_D = max(max(Dtemp(lines/3:2*lines/3,cols/3:2*cols/3)))+0.1;

    [EL,~,TotalPerf] = getBestFitEllipses(IClustNew,EL,NUMEllipses,area);
    
    if changes/area < 0.005 || ite > 40,
        %disp(sprintf('changes = %d ite = %d',changes,ite));
        break;
     end
    %[ok] = drawEllClusteting(IClustNew,EL,0,0);

end
TotalPerf;







