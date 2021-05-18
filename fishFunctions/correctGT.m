function [ GT ] = correctGT( GT)



GT(GT >= 1) = 1;
GT(GT < 1) = 0;


end

