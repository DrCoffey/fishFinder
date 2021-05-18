%Initilization of Ellipses
function [EL,NUMEllipses] = initEll(Cent,Rad,ELLSET)

NUMEllipses = length(ELLSET);
EL = [];
for val=1:length(ELLSET),
    id = ELLSET(val);
    EL(val).a = Rad(id);
    EL(val).b = Rad(id);
    if val == 1,
        EL(val).a = Rad(id);
        EL(val).b = Rad(id);
    end
    EL(val).C(1) = Cent(id,2);
    EL(val).C(2) = Cent(id,1);
    EL(val).phi = 0;
    EL(val).InArea = 0;
    EL(val).outPixels = 0;
    EL(val).tomh_area = 0;
    EL(val).tomh_enwsh = 0;
    EL(val).Label = val;
    EL(val).ELLSET = id;
end

end
