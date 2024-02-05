clc;
clear all;
close all;

%
%  This routine sets up the eigenvalue problem curl curl u = lambda u
% with the lowest order approximation for the Hodge

L = pi;
N = 30;
M = 30;

dx = L/N;
dy = L/M;

A = sparse(zeros(2*N*M-N-M));

nou = (N-1)*M;  % "nou" denotes the number of unknown u's

for i=1:N-1
    for j=1:M
        row = i + (j-1)*(N-1);
        if j>1
            col = i+(j-2)*(N-1);
            A(row,col) = -1/(dy*dy);
            col = nou + i + (j-2)*N;
            A(row,col) = 1/(dx*dy);
            col = nou + i+1 + (j-2)*N;
            A(row,col) = -1/(dx*dy);
        end
        col = row;
        A(row,col) = 2/(dy*dy);
        if j<M
            col = i + j*(N-1);
            A(row,col) = -1/(dy*dy);
            col = nou + i + (j-1)*N;
            A(row,col) = -1/(dx*dy);
            col = nou + i + 1 + (j-1)*N;
            A(row,col) = 1/(dx*dy);
        end
    end
end

for i=1:N
    for j=1:M-1
        row = nou + i +(j-1)*N;
        if i>1
            col = i-1 + (j-1)*(N-1);
            A(row,col) = 1/(dx*dy);
            col = i-1 + (j-1+1)*(N-1);
            A(row,col) = -1/(dx*dy);
            col = nou + i - 1+(j-1)*N;
            A(row,col) = -1/(dx*dx);
        end
        col = row;
        A(row,col) = 2/(dy*dy);
        if i<N
            col = i + (j-1)*(N-1);
            A(row,col) = -1/(dx*dy);
            col = i +(j-1+1)*(N-1);
            A(row,col) = 1/(dx*dy);
            col = nou + i + 1 + (j-1)*N;
            A(row,col) = -1/(dx*dx);
        end
    end
end

% opts.issym = 1;
% lambda = eigs(A,2*N*M-N-M,'sa',opts);
lambda = sort(eig(full(A)));

teller = 1;
while abs(lambda(teller))<1e-6, teller = teller + 1; end

plot(lambda,'.')

figure(2)
plot(lambda(teller:teller+9),'*')
xlim([0,10])
grid on
