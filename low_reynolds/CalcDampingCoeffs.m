function [ aCoeff ] = CalcDampingCoeffs( L_sq, dY, length , BC)
%CALCCOEFFS Summary of this function goes here
%   Detailed explanation goes here

    dYnorth = dY(:,1);
    dYsouth = dY(:,2);
    
    aCoeff.north = zeros(length,1);
    aCoeff.south = zeros(length,1);

    for i = 2:length-1

        aCoeff.north(i) = (-L_sq(i))/(dYnorth(i)*(dYnorth(i) - dYsouth(i)));
        aCoeff.south(i) = (-L_sq(i))/(dYsouth(i)*(dYnorth(i) - dYsouth(i)));

    end
    
    % Not sure about this
    aCoeff.south(2) = aCoeff.south(2) * BC(1);
    aCoeff.north(end-1) = aCoeff.north(end-1) * BC(2);
    
    aCoeff.point = 1 - aCoeff.north - aCoeff.south; 


end

