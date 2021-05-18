%Returns the equal area BestFitEllipse 
function [EL,area,p] = getBestFitEllipse(I,EL,val)
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
EL(val).Label = val;
EL(val).ELLSET = EL(val).ELLSET;

p = p1(:,1)+lines*cols*p1(:,2);

