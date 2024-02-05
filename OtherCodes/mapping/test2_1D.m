clear all
clf%ose all
clc

N = 7;

X = linspace(-1,1,1000);
PhiCoeff = poly(2*rand(1,N)-1);
UCoeff = polyder(PhiCoeff);
Phi = zeros(1,1000);
for n=0:N
    Phi = Phi + PhiCoeff(N+1-n)*X.^n;
end

plot(X,Phi)

U = zeros(1,1000);
for n=0:N-1
    U = U + UCoeff(N-n)*X.^n;
end

hold on
plot(X,U,'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xigl,wgl] = GLLnodes(N);

[hglgl,dhdxglgl] = LagrangeVal(xigl,N,1);

% xgl = xigl;
% dxdxigl = ones(size(xigl));

xgl = 1/3*(3/2+xigl).^2-13/12;
dxdxigl = 1+2/3*xigl;

% xgl = 1/4*(1+xigl).^3-1;
% dxdxigl = 3/4*(1+xigl).^2;

phi = zeros(1,N+1);
for n=0:N
    phi = phi + PhiCoeff(N+1-n)*xgl.^n;
end


A = diag(wgl);
B = zeros(N+1);
for p=1:N+1
    for i=1:N+1
        B(p,i) = wgl(p)*dhdxglgl(i,p)/dxdxigl(p);
    end
end

D = inv(A)*B;

u = D*phi';

hh = LagrangeVal(X,N,1);

uu = u'*hh;

% XX = X;

XX = 1/3*(3/2+X).^2-13/12;

% XX = 1/4*(1+X).^3-1;

plot(XX,uu,'--g')

plot(xgl,phi,'xc')