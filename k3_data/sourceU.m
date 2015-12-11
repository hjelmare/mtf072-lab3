function [Sp,Su] = sourceU(deltaY,length)

    Sp = zeros(length,1);
    Su = ones(length,1) .* deltaY;

end