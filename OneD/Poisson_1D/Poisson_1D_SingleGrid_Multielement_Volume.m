clear all
close all
clc

NrElementRange = 2;
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

dxxdXiXi = Jac*ones(1,kk*M);

for N=NrCellRange
[xi,w] = GLLnodes(N);

Xi = zeros(1,N*M+1);
for m=1:M
    ind = N*(m-1)+(1:N+1);
    Xi(ind) = xi+2*(m-1);
end

x = (Xi+1)*Jac;

dxdXi = Jac*ones(1,N*M+1);

phi = exp(x).*sin(pi*x);

f = exp(x).*(sin(pi*x)+pi*cos(pi*x));

%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

D = full(topology1D(N*M));

M0 = zeros(1,N*M+1);
for m=1:M
    ind = N*(m-1)+(1:N+1);
    M0(ind) = M0(ind) + w.*dxdXi(ind);
end
M0 = diag(M0);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(w.*e(i,:).*e(j,:)./dxdXi(1:N+1));
    end
end

M1 = kron(eye(M),M1); % Only possible for uniform grid !!!

L = -M1*D/M0*D'*M1;


% p_bc = [ -(-pi*exp(-1)) ; zeros(size(L,1)-2,1) ; pi*exp(1) ];
F = M1*diff(f)';% + p_bc;
    
%% boundary conditions
% ????????????????????????????

% % solve
% a = L\F;
% 
% a = a';

a = [ M0 D'*M1 ; M1*D zeros(N*M) ]\[zeros(N*M+1,1) ; F];

a = a(N*M+1+(1:N*M));

%%

[hh,ee]=MimeticpolyVal(xixi,N,1);
aa = zeros(size(xx));
for m=1:M
    ind1 = kk*(m-1)+(1:kk);
    ind2 = N*(m-1)+(1:N);
    aa(ind1) = a(ind2)'*ee/Jac;
end


clf
plot(xx,aa_ex,'r')
hold on
plot(xx,aa,'--b')
for i=1:N*M
    plot([x(i) x(i+1)],[a(i) a(i)]/(x(i+1)-x(i)),'g')%*dxdXi(i)
end
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
break
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