function [ GT2 ] = convertGTBoundarytoRegions( GT )

lines = size(GT,1);
cols = size(GT,2);
GT0 = GT;

for i=1:cols,
    GT0(1,i) = 1;
end
for i=1:lines,
    GT0(i,cols) = 1;
end
for i=1:cols,
    GT0(lines,i) = 1;
end
GTH0 = imfill(GT0,'holes');

GT1 = GT;
for i=1:cols,
    GT1(1,i) = 1;
end
for i=1:lines,
    GT1(i,1) = 1;
end
for i=1:cols,
    GT1(lines,i) = 1;
end
GTH1 = imfill(GT1,'holes');
GT2 = max(GTH0,GTH1);
GT2(GT == 1) = 0;

for i=1:cols,
    if GT2(2,i) == 0
        GT2(1,i) = 0;
    end
end
for i=1:lines,
    if GT2(i,cols-1) == 0
        GT2(i,cols) = 0;
    end
end
for i=1:lines,
    if GT2(i,2) == 0
        GT2(i,1) = 0;
    end
end
for i=1:cols,
    if GT2(lines-1,i) == 0
        GT2(lines,i) = 0;
    end
end



end

