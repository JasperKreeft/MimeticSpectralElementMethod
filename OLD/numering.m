clear; clc;

numRows    = 2;
numColumns = 3;

N = 3;


localnr = reshape(1:N*N,N,N);


globalnr = zeros(numColumns*N,numRows*N);
for j=1:numRows
    for i=1:numColumns
        globalnr((i-1)*N+(1:N),(j-1)*N+(1:N)) = ((i-1)+(j-1)*numColumns)*N*N+localnr;
    end
end

disp('globalnr = '); disp('  '); disp(flipud(globalnr'));


