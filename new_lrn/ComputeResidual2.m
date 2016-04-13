function res = ComputeResidual2(U,UCoeff,uSu,R,RCoeff,RSu, f, fCoeff, ...
    fSu,nj)

   %Extracting variables from struct variables
   aPU = UCoeff.point;
   aNU = UCoeff.north;
   aSU = UCoeff.south;

   aPR = RCoeff.point;
   aNR = RCoeff.north;
   aSR = RCoeff.south;
   
   aPf = fCoeff.point;
   aNf = fCoeff.north;
   aSf = fCoeff.south;

        
   UR = sum(abs(aPU(2:nj-1).*U(2:nj-1) - aSU(2:nj-1).*U(1:nj-2) - ...
       aNU(2:nj-1).*U(3:nj) - uSu(2:nj-1)));
   RR = sum(abs(aPR(2:nj-1).*R(2:nj-1) - aSR(2:nj-1).*R(1:nj-2) - ...
       aNR(2:nj-1).*R(3:nj) - RSu(2:nj-1)));
   
   fR = sum(abs(aPf(2:nj-1).*f(2:nj-1) - aSf(2:nj-1).*f(1:nj-2) - ...
       aNf(2:nj-1).*f(3:nj) - fSu(2:nj-1)));
   
    
    res = UR + RR + fR; 


end