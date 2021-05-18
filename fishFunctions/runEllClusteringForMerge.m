%Performs clustering based on ellipse set EL 
function [EL,IClustNew,Dtemp,TotalPerf] = runEllClusteringForMerge(EL,ELLSET,IClust,area)
INF = 100000000000000000;
IClustNew = IClust;
lines = size(IClust,1);
cols = size(IClust,2);
NUMEllipses = length(ELLSET);
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
                d = INF*ones(1,max(ELLSET));
                for kid=1:NUMEllipses,
                    k = ELLSET(kid);
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
                IClustNew(i,j) = pos; %difference
            end
        end
    end
    Thresh_D = max(max(Dtemp(lines/3:2*lines/3,cols/3:2*cols/3)))+0.1;

    [EL,~,TotalPerf] = getBestFitEllipsesForMerge(IClustNew,EL,ELLSET,area);
    if ite > 1 || changes == 0,
        if changes/area < 0.005 || ite > 40,
            %disp(sprintf('changes = %d ite = %d',changes,ite));
            break;
        end
    end
    %[ok] = drawEllClusteting(IClustNew,EL,0,0);

end
TotalPerf;

%[ok] = drawEllClusteting(IClustNew(lines/3+1:2*lines/3,cols/3+1:2*cols/3),EL,lines/3,cols/3);

%[ok] = drawDistEllClusteting(Dtemp(lines/3+1:2*lines/3,cols/3+1:2*cols/3),EL,lines/3,cols/3);



function [EL,area,TotalPerf] = getBestFitEllipsesForMerge(I,EL,ELLSET,area)
p = [];
%s = 0;
for kid=1:length(ELLSET),
    k = ELLSET(kid);
    [EL,~,p1] = getBestFitEllipseForMerge(I,EL,k);
    p = union(p,p1);
   % s = s+length(p1);
end
%TotalPerf = s/area;
TotalPerf = size(p,1)/area;






%Returns the equal area BestFitEllipse 
function [EL,area,p] = getBestFitEllipseForMerge(I,EL,val)
%I = imrotate(I,30,'nearest','loose');
BW = I == val;
BW0 = I > 0;    
lines = size(I,1);
cols = size(I,2);

[x y] = meshgrid(1:max(lines,cols),1:max(lines,cols));
X0 = EL(val).C(1);
Y0 = EL(val).C(2);
el=((x-X0)/EL(val).a).^2+((y-Y0)/EL(val).b).^2<=1;
% radE = sqrt(EL(val).a*EL(val).b);
% se = strel('disk',ceil(0.05*radE));
% el = imerode(el,se);

el = rotateAround(el,Y0,X0,EL(val).phi,'nearest');
el = el(1:lines,1:cols);
el = min(el,BW0);
BW1 = max(BW,el);


stats = regionprops(double(BW1), 'Area','Centroid','MajorAxisLength','MinorAxisLength','Orientation');
if sum(sum(BW)) == 0,
    p = [];
    area = 0;
    return;
end
area = stats.Area;
%RealArea = sum(BW(:));
%area = mean([RealArea area]);

C = stats.Centroid;
e =  stats.MajorAxisLength / stats.MinorAxisLength;
X0 = C(1);
Y0 = C(2);
phi = stats.Orientation;
%pi a b = area
% a/b = e
% a = e*b
% a^2 = e*area/pi

a = sqrt(e*area/pi);
b = a/e;
% 
% apoX = round(min(1,X0-2*a-2));
% eosX = round(max(cols,X0+2*a+2));
% 
% apoY = round(min(1,X0-2*b-2));
% eosY = round(max(lines,X0+2*b+2));
% 
% [x y] = meshgrid(apoX:eosX,apoY:eosY);

%[x y] = meshgrid(1:max(lines,cols),1:max(lines,cols));
%
el=((x-X0)/a).^2+((y-Y0)/b).^2<=1;

el = rotateAround(el,Y0,X0,phi,'nearest');
el = el(1:lines,1:cols);

p1 = [];
p2 = [];
[p1(:,1) p1(:,2)] = find(el == 1 & BW == 1);
[p2(:,1) ~] = find(el == 1 | BW == 1);

tomh_area = size(p1,1) / area;
tomh_enwsh = size(p1,1) / size(p2,1);

EL(val).a = a;
EL(val).b = b;
EL(val).C = C;
EL(val).phi = phi;
EL(val).InArea = size(p1,1);
EL(val).outPixels = size(p2,1) - size(p1,1);
EL(val).tomh_area = tomh_area;
EL(val).tomh_enwsh = tomh_enwsh;
EL(val).Label = val; %difference
EL(val).ELLSET = EL(val).ELLSET;

p = p1(:,1)+lines*cols*p1(:,2);



