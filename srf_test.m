% clear workspace; clc;

% the script do a nummerical test of srf algorithm
% generate the sparase matrix
% set n1 = n2
n1 = 50; n2 = 50; r = 5;
ML = rand(n1,r);    MR = rand(r,n2);
M = ML*MR;

% degree of freedom and number of enteries
dof = r*(n1+n2-r);   m = 2*dof;
% enteries postitions
s = rand(1,n1*n2);
[~,I] = sort(s,2);
omega = zeros(size(M)); Me = zeros(size(M));
row = ceil(I/n2); col = rem(I,n1)+1;
for oi = 1:m
    omega(row(oi),col(oi)) = 1;
    Me(row(oi),col(oi)) = M(row(oi),col(oi));
end

Xp = srf(Me,omega);
err = norm(M-Xp,'fro')/norm(M,'fro') ;