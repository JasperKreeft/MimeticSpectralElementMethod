clear all
close all
clc

NrElementRange = 2:6;
NrCellRange = 2;

kk = 1000;
xixi=linspace(-1,1,kk);


L2 = zeros(1,size(NrCellRange,2)*size(NrElementRange,2));
C  = zeros(1,size(NrCellRange,2)*size(NrElementRange,2));
n=0;

for M=NrElementRange

Jac = 1/(2*M);
    
XiXi = zeros(1,kk*M);
for m=1:M
    ind = kk*(m-1)+(1:kk);
    XiXi(ind) = xixi+2*(m-1);
end
xx = (XiXi+1)*Jac;
aa_ex = exp(xx).*sin(pi*xx);

for N=NrCellRange
[xi,w] = GLLnodes(N);

Xi = zeros(1,N*M+1);
for m=1:M
    ind = N*(m-1)+(1:N+1);
    Xi(ind) = xi+2*(m-1);
end

x = (Xi+1)*Jac;

phi = exp(x).*sin(pi*x);

f = exp(x).*((1-pi^2)*sin(pi*x)+2*pi*cos(pi*x));

[h,e] = MimeticpolyVal(xi,N,1);

G = full(topology1D(N*M));

A1 = zeros(N);
for i=1:N
    for j=1:N
        A1(i,j) = sum(w.*e(i,:).*e(j,:))/Jac;
    end
end

A = kron(eye(M),A1);

Li = -G'*A*G;

% B = zeros(N*M+1,N*M);
% for i=1:M
%     ind = N*(i-1)+(1:N);
%     B(N*(i-1)+1,ind) = -e(:,1);
%     B(N*i+1,ind)     =  e(:,N+1);
% end
% 
% Lb = B*G;

% L = -Li+Lb;

% L = L*((2*M)^2); % Jacobiaan

L = Li;

A0 = zeros(1,N*M+1);
for m=1:M
    ind = N*(m-1)+(1:N+1);
    A0(ind) = A0(ind) + w*Jac;
end
A0 = diag(A0);

F = A0*f';

%% boundary conditions
L_full = L;

L(N*M+1,:) = [];
L(1,:)     = [];
L(:,N*M+1) = [];
L(:,1)     = [];

F_full = F;
F = F(2:N*M);

% solve
a = L\F;

a = [ 0  a' 0 ];

%%

[hh,ee]=MimeticpolyVal(xixi,N,1);
aa = zeros(size(xx));
for m=1:M
    ind1 = kk*(m-1)+(1:kk);
    ind2 = N*(m-1)+(1:N+1);
    aa(ind1) = a(ind2)*hh;
end


clf
plot(xx,aa_ex,'r')
hold on
plot(xx,aa,'--b')
plot(x,a,'x')
grid

% plot(x,phi)
% hold on
% plot(x,a,'r')

%%
n=n+1
L2(n) = sqrt( 2/kk*sum((aa-aa_ex).^2) );
C(n) = cond(L);
end
end
nhf = floor(n/2);

figure
% semilogy(NrCellRange,L2,'-')
loglog(NrElementRange,L2,'-')
title('L_2 convergence')
grid on

figure
% Line2 = 0.05*NrCellRange.^3;
% loglog(NrCellRange,Line2,'--g')
% str2 = 'y=0.05x^3';
% hold on
% rate = (log(C(n))-log(C(nhf)))/(log(NrCellRange(n))-log(NrCellRange(nhf)));
% Line = C(n)/NrCellRange(n)^rate*NrCellRange.^rate;
% loglog(NrCellRange,Line,'r')
% hold on
% str = ['y=' num2str(C(n)/NrCellRange(n)^rate,3) 'x^{' num2str(rate,3) '}'];
% loglog(NrCellRange,C,'-');
loglog(NrElementRange,C,'-');
grid on
y=ylim;
% xlim([NrCellRange(1) NrCellRange(n)])
ylim([1 y(2)])
% legend(str2,'Condition number',4)
% legend(str2,'Condition number',str,4)