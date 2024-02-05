% Dirichlet
% 2 5 5 8 10 10 13 13 17 17 18 20 20 25 25 26 26 29 29 32 34 34 37 37 40...
for m=1:10
    for n=1:10
        ind = n+(m-1)*10;
        EexD(ind) = n^2+m^2;
    end
end
EexD = sort(EexD)';


% Neumann
% 1 1 2 4 4 5 5 8 9 9 10 10 13 13 16 16 17 17 18 20 20 25 25 25 25 26 26...
for m=0:80
    for n=0:80
        ind = (n+1)+m*81;
        EexN(ind) = n^2+m^2;
    end
end
EexN = sort(EexN)';
EexN(1) = [];