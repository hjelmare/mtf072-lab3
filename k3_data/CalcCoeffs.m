function [ aCoeff ] = CalcCoeffs( sigma, dY, viscosity, vist, length )
%CALCCOEFFS Summary of this function goes here
%   Detailed explanation goes here

aCoeff.point = zeros(length,1);
aCoeff.north = zeros(length,1);
aCoeff.south = zeros(length,1);

for i = 1:length
    aCoeff.point(i) = (((viscosity + vist(i+1)/sigma)/dY(i+1)) + ((viscosity + vist(i-1)/sigma)/dY(i-1)));
    aCoeff.north(i) = ((viscosity + vist(i+1)/sigma)/dY(i+1));
    aCoeff.south(i) = ((viscosity + vist(i-1)/sigma)/dY(i-1));
end



end

