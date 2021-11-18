function [Entropy] = EntropyCalc(Type,sourceVar)
%
if strcmp(Type, 'Gaussian')
    Entropy = 0.5 * log(2*pi*exp(1) * sourceVar); 
elseif  strcmp(Type, 'Exp')
    lambda = 1/sqrt(sourceVar);
    Entropy = 1 - log(lambda); 
elseif  strcmp(Type, 'Laplace')
    b = sqrt(sourceVar/2);
    Entropy = 1 + log(2*b);    
elseif strcmp(Type, 'Unifrom')
    Delta = sqrt(3*sourceVar);
    Entropy = log(Delta); 
end
end

