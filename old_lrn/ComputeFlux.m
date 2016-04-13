function F = ComputeFlux(U,deltaY,nj)

    F = sum(U(2:nj-1) .* deltaY(2:nj-1));


end