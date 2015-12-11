function [ aCoeff ] = CalcCoeffs( sigma, dY, viscosity, vist, length )
%CALCCOEFFS Summary of this function goes here
%   Detailed explanation goes here

dYnorth = dY(:,1);
dYsouth = dY(:,2);

aCoeff.point = zeros(length,1);
aCoeff.north = zeros(length,1);
aCoeff.south = zeros(length,1);

for i = 1:length
    aCoeff.point(i) = (((viscosity + vist(i+1)/sigma)/dYnorth(i)) + ((viscosity + vist(i-1)/sigma)/dYsouth(i)));
    aCoeff.north(i) = ((viscosity + vist(i+1)/sigma)/dYnorth(i));
    aCoeff.south(i) = ((viscosity + vist(i-1)/sigma)/dYsouth(i));
end



end

