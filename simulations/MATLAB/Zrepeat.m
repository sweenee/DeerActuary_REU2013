fp=fopen('Z.csv','w');
fprintf(fp,'Z\n');

for i=1:1000
    [Z, X] = DeerPop(10, 300, 30, 20, 3000000, 1);
     fprintf(fp,'%f\n',Z(300+1)); %records last position Z
end

fclose(fp); 