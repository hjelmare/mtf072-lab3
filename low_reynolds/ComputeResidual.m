function R = ComputeResidual(U,UCoeff,uSu,k,kCoeff,kSu,eps,epsCoeff, ...
    epsSu,nj)

   %Extracting variables from struct variables
   aPU = UCoeff.point;
   aNU = UCoeff.north;
   aSU = UCoeff.south;

   aPk = kCoeff.point;
   aNk = kCoeff.north;
   aSk = kCoeff.south;

   aPeps = epsCoeff.point;
   aNeps = epsCoeff.north;
   aSeps = epsCoeff.south;

        
   UR = sum(abs(aPU(2:nj-1).*U(2:nj-1) - aSU(2:nj-1).*U(1:nj-2) - ...
       aNU(2:nj-1).*U(3:nj) - uSu(2:nj-1)));
   kR = sum(abs(aPk(2:nj-1).*k(2:nj-1) - aSk(2:nj-1).*k(1:nj-2) - ...
       aNk(2:nj-1).*k(3:nj) - kSu(2:nj-1)));
   epsR = sum(abs(aPeps(2:nj-1).*eps(2:nj-1) - aSeps(2:nj-1).*eps(1:nj-2)-...
       aNeps(2:nj-1).*eps(3:nj) - epsSu(2:nj-1)));        
   
    
    R = UR + kR + epsR; 


end