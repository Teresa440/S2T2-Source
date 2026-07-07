clear all; clc;

g = @(x) 1+9.*sum(x(:,2:end),2)./(size(x,2)-1);
MultiObj.fun = @(x) [x(:,1), g(x).*(1-sqrt(x(:,1)./g(x)))];
MultiObj.nVar = 30; 
MultiObj.var_min = zeros(1,MultiObj.nVar);
MultiObj.var_max = ones(1,MultiObj.nVar);

RVEA(MultiObj.fun, MultiObj.var_min, MultiObj.var_max);
