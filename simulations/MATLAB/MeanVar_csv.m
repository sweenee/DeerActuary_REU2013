function [  ] = MeanVar_csv( alpha, alpha0, dalpha, gamma, gamma0, dgamma, P, P0, dP )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

for i=alpha0:dalpha:alpha
    for j=gamma0:dgamma:gamma
        for k=P0:dP:P
            fp=fopen('meanvar.csv','a');
            fprintf(fp,'Xmean\n', fp,'Xvar\n', fp,'Mmean\n', fp,'Mvar\n');
            [Xmean, Xvar, Mmean, Mvar] = MeanVar_DeerInsSim(i, j, k);
            fprintf(fp,'%f\n',Xmean,fp,'%f\n',Xvar,fp,'%f\n',Mmean,fp,'%f\n',Mvar);
            fclose(fp);
        end
    end
end