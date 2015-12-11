function T = GaussSeidel(T,aCoeff)

    rows = length(T);
    
    %Defining source term Su
    Su = zeros(rows,1); %No source
    
    %Extracting variables from struct variable
    aP = aCoeff.point;
    aN = aCoeff.north;
    aS = aCoeff.south;

    
    %Performing Gauss-Seidel update
      for i = 2:rows-1
          %Calculating parts of update equation
          Tn = T(i+1);
          Ts = T(i-1);
            
          %Performing update on T(i,j)
          T(i) = (aN(i)*Tn + aS(i)*Ts + Su(i))/aP(i);
                  
      end
    
end