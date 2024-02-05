clear all
clf% close all
% clc

NrCellRange = 6;
% c = 0;

L2 = zeros(size(NrCellRange));
C  = zeros(size(NrCellRange));
n=0;
for N=NrCellRange
    
[xgl,wgl] = GLLnodes(N);

X = xgl + 0.2*sin(pi*xgl);

Jac = 1+0.2*pi*cos(pi*xgl);

[h,e]=MimeticpolyVal(xgl,N,1);

NG = topology1D(N);

M0 = diag(wgl.*Jac);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(wgl.*e(i,:).*e(j,:)./Jac);
    end
end

% f = exp(X).*((1-pi^2)*sin(pi*X)+2*pi*cos(pi*X));
f = -pi^2*sin(pi*X);

L = -NG'*M1*NG;
F = M0*f';

%% boundary conditions
L_full = L;

L(N+1,:) = [];
L(1,:)   = [];
L(:,N+1) = [];
L(:,1)   = [];

F_full = F;
F = F(2:N);

% solve
a = L\F;

a = [ 0  a' 0 ];

%%

kk = 1000;
xx=linspace(-1,1,kk);
yy = xx + 0.2*sin(pi*xx);
[hh,ee]=MimeticpolyVal(xx,N,1);

aa    = a*hh;
% aa_ex = exp(yy).*sin(pi*yy);
aa_ex = sin(pi*yy);

clf
plot(yy,aa_ex,'r')
hold on
plot(yy,aa,'--b')
plot(X,a,'x')
grid

%%
n=n+1;
L2(n) = sqrt( 2/kk*sum((aa-aa_ex).^2) );
C(n) = cond(L);

end
nhf = floor(n/2);

% figure
% semilogy(NrCellRange,L2,'-')
% title('L_2 convergence')
% grid on

% figure
% Line2 = 0.05*NrCellRange.^3;
% loglog(NrCellRange,Line2,'--g')
% str2 = 'y=0.05x^3';
% hold on
% rate = (log(C(n))-log(C(nhf)))/(log(NrCellRange(n))-log(NrCellRange(nhf)));
% Line = C(n)/NrCellRange(n)^rate*NrCellRange.^rate;
% % loglog(NrCellRange,Line,'r')
% str = ['y=' num2str(C(n)/NrCellRange(n)^rate,3) 'x^{' num2str(rate,3) '}'];
% loglog(NrCellRange,C,'-');
% grid on
% y=ylim;
% xlim([NrCellRange(1) NrCellRange(n)])
% ylim([1 y(2)])
% legend(str2,'Condition number',4)
% % legend(str2,'Condition number',str,4)