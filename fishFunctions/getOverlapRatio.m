%It gets the overlap of a set of Ellipses
function [ overlap,overlapMax ] = getOverlapRatio(EL,SET,I)

lines = size(I,1);
cols = size(I,2);
overlap = zeros(1,length(SET));

[x y] = meshgrid(1:max(lines,cols),1:max(lines,cols));
Map = zeros(lines,cols);
for i=1:length(SET),
    val = SET(i);
    X0 = EL(val).C(1);
    Y0 = EL(val).C(2);
    el=((x-X0)/EL(val).a).^2+((y-Y0)/EL(val).b).^2<=1;
    % radE = sqrt(EL(val).a*EL(val).b);
    % se = strel('disk',ceil(0.05*radE));
    % el = imerode(el,se);
    el = rotateAround(el,Y0,X0,EL(val).phi,'nearest');
    el = el(1:lines,1:cols);
    Map = Map+el;
end

for i=1:length(SET),
    val = SET(i);
    X0 = EL(val).C(1);
    Y0 = EL(val).C(2);
    el=((x-X0)/EL(val).a).^2+((y-Y0)/EL(val).b).^2<=1;
    % radE = sqrt(EL(val).a*EL(val).b);
    % se = strel('disk',ceil(0.05*radE));
    % el = imerode(el,se);
    el = rotateAround(el,Y0,X0,EL(val).phi,'nearest');
    el = el(1:lines,1:cols);
    E = sum(el(:));
    [a,b] = find(Map == 1 & el == 1);
    overlap(i) =  (E-length(a)) / E;
end

overlapMax = max(overlap);%(totalArea-sum(objectArea)) / totalArea;
overlap = mean(overlap.^2);%(totalArea-sum(objectArea)) / totalArea;

if overlapMax < 0,
    overlap = 0;
end
%overlap = 0;
end

%It gets the overlap of a set of Ellipses
function [ overlap,totalArea ] = getOverlapRatio2(EL,SET)

EllArea = zeros(1,length(SET));
objectArea = zeros(1,length(SET));
overlap = zeros(1,length(SET));
for i=1:length(SET),
    k = SET(i);
    objectArea(i) = EL(k).InArea;
    EllArea(i) = pi*EL(k).a*EL(k).b;
    overlap(i) = (EllArea(i)-objectArea(i))/ EllArea(i);
end

totalArea = sum(EllArea);
overlap = max(overlap);%(totalArea-sum(objectArea)) / totalArea;

if overlap < 0,
    overlap = 0;
end
%overlap = 0;
end





