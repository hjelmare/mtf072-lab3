function [ aCoeff ] = CalcCoeffs( sigma, dY, deltaY, viscosity, vist, ...
    Sp, length , BC)
%CALCCOEFFS Summary of this function goes here
%   Detailed explanation goes here

    dYnorth = dY(:,1);
    dYsouth = dY(:,2);
    
    aCoeff.north = zeros(length,1);
    aCoeff.south = zeros(length,1);

    for i = 2:length-1
        
        vistNorth = (vist(i+1) - vist(i))*deltaY(i) / (2*dYnorth(i));
        vistSouth = (vist(i) - vist(i-1))*deltaY(i) / (2*dYsouth(i));
        
        aCoeff.north(i) = ((viscosity + vistNorth/sigma)/dYnorth(i));
        aCoeff.south(i) = ((viscosity + vistSouth/sigma)/dYsouth(i));

    end
    
    aCoeff.south(2) = aCoeff.south(2) * BC(1);
    aCoeff.north(end-1) = aCoeff.north(end-1) * BC(2);
    
    aCoeff.point = aCoeff.north + aCoeff.south - Sp; 


end

