function [Sp,Su] = sourceK(deltaY,dudy,vist,epsi,k,length)

    Su = (vist .* (dudy).^2) .* deltaY;
    
    Sp = (-epsi./ k) .* deltaY;
    
    

end