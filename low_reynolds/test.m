clear all
clc

T  = [10 15 20 25 30 35 20 20 20 20 20 20 10];

aCoeff.point = ones(length(T),1);
aCoeff.north = ones(length(T),1) * 0.5;
aCoeff.north(end-1) = 0;
aCoeff.south = aCoeff.north(end:-1:1);

T_store = [];

for i = 1:40
    T_new = GaussSeidel(T,zeros(size(T)),aCoeff);
    T_store = [T_store ; T_new];
    T = T_new;
end

contourf(T_store)