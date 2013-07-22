function [ Xtrue, Mem ] = DeerInsSim(  T, N, r1, h, F, alpha, rho, beta, P, gamma, g )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dt = T/N;
s = (0:dt:T);

clf;

%Deer Population
subplot(1, 2, 1);
[Ztrue, Xtrue, Xem] = DeerPop_Beth( T, N, r1, h, F, alpha );
plot(s, Xtrue, 'g-')
hold on
plot(s, Xem, 'r*')
%plot(s, Xmil, 'b*')
title('Deer Population');
legend('Xtrue', 'Xem')

%Insurance Payout
subplot(1, 2, 2);
[Mem] = InsPayout(T, N, r1, h, F, alpha,  rho, beta, P, gamma, g);
plot(s, Mem,'r*')
%hold on
%plot( s, Mmil, 'b*')
title('Insurance Payout');
legend( 'Mem')


