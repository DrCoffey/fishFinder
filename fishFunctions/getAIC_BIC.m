% AIC - BIC computation 
 

function [AIC,BIC,RES,bestAICBIC,SI] = getAIC_BIC(nCompl,TotalPerf,NUMEllipses,AICBIC_SELECTION,IClust,EL,ELLSET)
CONST = 1;
MODEL_PAR = 1;
%[ overlap ] = getOverlapRatio(EL,[1:NUMEllipses]);
SI = [0 0 0];
%[ overlap,~ ] = getOverlapRatio(EL,ELLSET,IClust);

AIC = nCompl*log((1-TotalPerf))+2*CONST*MODEL_PAR*NUMEllipses;
BIC = nCompl*log((1-TotalPerf))+MODEL_PAR*CONST*NUMEllipses*log(nCompl);

RES = BIC;
if AICBIC_SELECTION == 1,
    RES = AIC;
end

if TotalPerf > 0.97
    maxAIC = nCompl*log(1-0.9999)+2*CONST*MODEL_PAR*NUMEllipses;
    maxBIC = nCompl*log(1-0.9999)+MODEL_PAR*CONST*NUMEllipses*log(nCompl);    
elseif TotalPerf > 0.93
    maxAIC = nCompl*log(1-0.98)+2*CONST*MODEL_PAR*NUMEllipses;
    maxBIC = nCompl*log(1-0.98)+MODEL_PAR*CONST*NUMEllipses*log(nCompl);
else
    maxAIC = nCompl*log(1-0.96)+2*CONST*MODEL_PAR*NUMEllipses;
    maxBIC = nCompl*log(1-0.96)+MODEL_PAR*CONST*NUMEllipses*log(nCompl);
end

bestAICBIC = maxBIC;
if AICBIC_SELECTION == 1,
    bestAICBIC = maxAIC;
end












