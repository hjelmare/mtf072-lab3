function field = GaussSeidel(field, Su,aCoeff)

    rows = length(field);
    
    %Extracting variables from struct variable
    aP = aCoeff.point;
    aN = aCoeff.north;
    aS = aCoeff.south;

    
    %Performing Gauss-Seidel update
      for i = 2:rows-1
          %Calculating parts of update equation
          Tn = field(i+1);
          Ts = field(i-1);
            
          %Performing update on T(i,j)
          field(i) = (aN(i)*Tn + aS(i)*Ts + Su(i))/aP(i);
                  
      end
    
end