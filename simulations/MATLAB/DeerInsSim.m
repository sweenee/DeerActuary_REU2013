function [  ] = DeerInsSim(  T, N, r1, h, F, alpha )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dt = T/N;
s = (0:dt:T);

clf;

subplot(1, 2, 1);
[Ztrue, Xtrue, Xem, Xmil] =DeerPop_Beth( T, N, r1, h, F, alpha );
plot(s, Xtrue, 'g-')
hold on
plot(s, Xem, 'r*')
plot(s, Xmil, 'b*')
title('Deer Population');
legend('Xtrue', 'Xem', 'Xmil')

subplot(1, 2, 2);
[Mem, Mmil] = InsPayout(T, N, r1, h, F, alpha);
plot(s, Mem,'r*',  s, Mmil, 'b*')
title('Insurance Payout');
legend( 'Mem', 'Mmil')


