function [ Xmean, Xvar, Mmean, Mvar ] = MeanVar_DeerInsSim( alpha, gamma, P )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

T = 10;
N = 25000;
r1 = ;
h = ;
F = ;
rho = ;
beta = ;
g = ;

for i=1:1000
    [Xtrue, Mem] = DeerInsSim(  T, N, r1, h, F, alpha, rho, beta, P, gamma, g );
end

Xmean = mean(Xtrue);
Xvar = var(Xtrue);
Mmean = mean(Mem);
Mvar = var(Mem);


