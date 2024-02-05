clear all
close all
clc

NrCellRange = 4;

kk = 1000;
xixi=linspace(-1,1,kk);

xx = (xixi+1)/2;
% uu_ex = exp(xx).*sin(pi*xx);
uu_ex = exp(xx).*cos(pi*xx);

dxxdxixi = 1/2*ones(1,kk);

for N=NrCellRange
[xi,w] = GLLnodes(N);

dxdxi = 1/2;

x = (xi+1)*dxdxi;

% phi = exp(x).*sin(pi*x);
% f = exp(x).*(sin(pi*x)+pi*cos(pi*x));
phi = exp(x).*cos(pi*x);
f = exp(x).*(cos(pi*x)-pi*sin(pi*x));

[h,e] = MimeticpolyVal(xi,N,1);

D = full(topology1D(N));


M0 = w*dxdxi;
M0 = diag(M0);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(w.*e(i,:).*e(j,:)/dxdxi);
    end
end

L = [ M0     D'*M1
      M1*D  zeros(N) ];

F = [ zeros(N+1,1)
       M1*diff(f)' ];

%% left flux boundary condition
% qL = pi; ubR = 0;
qL = 1; ubR = -exp(1);

F = F - L(:,1)*qL;
F = F + [zeros(N,1) ; ubR ; zeros(N,1)];

L(:,1) = [];
L(1,:) = [];

F(1) = [];

qu = L\F;

q = [qL ; qu(1:N)];
u = qu(N+1:2*N)';

%% right flux boundary condition
% % qR = -exp(1)*pi;
% qR = -exp(1); ubL = 1;
% 
% F = F - L(:,N+1)*qR;
% F = F + [ -ubL ; zeros(2*N,1) ];
% 
% L(:,N+1) = [];
% L(N+1,:) = [];
% 
% F(N+1) = [];
% 
% qu = L\F;
% 
% q = [qu(1:N) ; qR];
% u = qu(N+1:2*N)';

%%

[hh,ee]=MimeticpolyVal(xixi,N,1);
uu(1:kk) = u(1:N)*ee/dxdxi;

clf
plot(xx,uu_ex,'r')
hold on
plot(xx,uu,'--b')
for i=1:N
    plot([x(i) x(i+1)],[u(i) u(i)]/(x(i+1)-x(i)),'g')
end
grid on
legend('exact','numerical','discrete solution')

end % for N