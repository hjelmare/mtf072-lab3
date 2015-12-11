function [T,epsilon] = GaussSeidel(T,x,y,deltaX,deltaY,T1,c1,c2,kFactor)

    [rows,cols] = size(T);  
    
    %Predefining coefficients that handle boundary conditions
    eCoeff = ones(rows,cols);
    eCoeff(:,cols-1) = 2;
    wCoeff = ones(rows,cols);
    wCoeff(:,2) = 2;
    nCoeff = ones(rows,cols);
    nCoeff(rows-1,:) = 0;   % this is originally 0
    sCoeff = ones(rows,cols);
    sCoeff(2,:) = 2;
    
    %Performing Gauss-Seidel update
    maxR = 0;
    for j = 2:cols-1
        for i = 2:rows-1
            %Calculating parts of update equation
            dXeast = (x(j+1) - x(j));
            dXwest = (x(j) - x(j-1));
            dYnorth = (y(i+1) - y(i));
            dYsouth = (y(i) - y(i-1));
            Te = T(i,j+1);
            Tw = T(i,j-1);
            Tn = T(i+1,j);
            Ts = T(i-1,j);
            Tp = T(i,j);
            ae = (deltaY(i)/dXeast) * 2 *(1 + 20*Te/T1)*kFactor * eCoeff(i,j);
            aw = (deltaY(i)/dXwest) * 2 *(1 + 20*Tw/T1)*kFactor * wCoeff(i,j);
            an = (deltaX(j)/dYnorth) * 2 *(1 + 20*Tn/T1)*kFactor * nCoeff(i,j);
            as = (deltaX(j)/dYsouth) * 2 *(1 + 20*Ts/T1)*kFactor * sCoeff(i,j);
            bp = -15*c2*Tp*deltaX(j)*deltaY(i);
            bu = 15*c1*deltaX(j)*deltaY(i);
            ap = ae + aw + an + as - bp;
            
            %Performing update on T(i,j)
            T(i,j) = (ae*Te + aw*Tw + an*Tn + as*Ts + bu)/ap;
            
            %Calculating maximum residual
            R = abs(ap*Tp - ae*Te - aw*Tw - an*Tn - as*Ts - bu);
            if (R > maxR)
                maxR = R;
            end            
        end
    end
    
    epsilon = maxR;

end