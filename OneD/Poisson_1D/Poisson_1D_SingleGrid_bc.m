clear all
close all
clc

NrCellRange = 4:14;

L2 = zeros(size(NrCellRange));
C  = zeros(size(NrCellRange));
n=0;
for N=NrCellRange

[xgl,wgl] = GLLnodes(N);

Jac = 1/2;

[h,e]=MimeticpolyVal(xgl,N,1);

NG = topology1D(N);

M0 = diag(wgl*Jac);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(wgl.*e(i,:).*e(j,:))/Jac;
    end
end

ygl = (1+xgl)/2;
% f = exp(ygl).*((1-pi^2)*sin(pi*ygl)+2*pi*cos(pi*ygl));
f = exp(ygl).*((1-pi^2)*cos(pi*ygl)-2*pi*sin(pi*ygl));

L = -NG'*M1*NG;
F = M0*f';

%% Dirichlet boundary conditions
% % aL = 0; aR = 0;
% aL = 1; aR = -exp(1);
% 
% L_full = L;
% 
% F_full = F - L(:,1)*aL - L(:,N+1)*aR;
% F = F_full(2:N);
% 
% L(N+1,:) = [];
% L(1,:)   = [];
% L(:,N+1) = [];
% L(:,1)   = [];
% 
% % solve
% a = L\F;
% 
% a = [ aL  a' aR ];

%% left Neumann boundary conditions
% % aL = 0; aR = 0;
% qL = 1; aR = -exp(1);
% 
% L_full = L;
% 
% F_full = F - L(:,N+1)*aR - [ -qL ; zeros(N,1) ];
% F = F_full(1:N);
% 
% L(N+1,:) = [];
% L(:,N+1) = [];
% 
% % solve
% a = L\F;
% 
% a = [ a' aR ];

%% right Neumann boundary conditions
% aL = 0; aR = 0;
aL = 1; qR = -exp(1);

L_full = L;

F_full = F - L(:,1)*aL - [ zeros(N,1) ; qR ];
F = F_full(2:N+1);

L(1,:) = [];
L(:,1) = [];

% solve
a = L\F;

a = [ aL a' ];

%%

kk = 1000;
xx=linspace(-1,1,kk);
yy=(1+xx)/2;
[hh,ee]=MimeticpolyVal(xx,N,1);

aa    = a*hh;
% aa_ex = exp(yy).*sin(pi*yy);
aa_ex = exp(yy).*cos(pi*yy);

clf
plot(yy,aa_ex,'r')
hold on
plot(yy,aa,'--b')
plot(ygl,a,'x')
grid

%%
n=n+1;
L2(n) = sqrt( 2/kk*sum((aa-aa_ex).^2) );
C(n) = cond(L);

end
nhf = floor(n/2);

figure
semilogy(NrCellRange,L2,'-')
title('L_2 convergence')
grid on

figure
Line2 = 0.05*NrCellRange.^3;
loglog(NrCellRange,Line2,'--g')
str2 = 'y=0.05x^3';
hold on
rate = (log(C(n))-log(C(nhf)))/(log(NrCellRange(n))-log(NrCellRange(nhf)));
Line = C(n)/NrCellRange(n)^rate*NrCellRange.^rate;
% loglog(NrCellRange,Line,'r')
str = ['y=' num2str(C(n)/NrCellRange(n)^rate,3) 'x^{' num2str(rate,3) '}'];
loglog(NrCellRange,C,'-');
grid on
y=ylim;
xlim([NrCellRange(1) NrCellRange(n)])
ylim([1 y(2)])
legend(str2,'Condition number',4)
% legend(str2,'Condition number',str,4)