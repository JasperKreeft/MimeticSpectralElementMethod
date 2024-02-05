%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Adaptive 1D Poisson problem                                             %
% solved with Mimetic Spectral Element Method                             %
% Case: dd*a=f, where a is a 1-form (volume form)                         %
%                                                                         %
% Written by Jasper Kreeft - 2011                                         %
% Contact: j.j.kreeft@tudelft.nl                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings

% X = [ 0 0.5 0.75 1 ];
% X = linspace(0,1,4);
% X = [ 0 0.25 0.5 0.625 0.75 0.875 1.0 ];
X = [ 0 0.625 0.75 0.875 1.0 ];
% X = [0 1];

H = length(X)-1; % number of spectral elements
N = 6;           % number of cells in an element

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation

[xi,w] = GLLnodes(N);

Xi       = zeros(1,N*H+1);
x        = zeros(1,N*H+1);
Jacobian = zeros(H,N+1);
for m=1:H
    ind = N*(m-1)+(1:N+1);
    Xi(ind) = xi+2*(m-1);

    Jacobian(m,:) = (1/2)*(X(m+1)-X(m))*ones(1,N+1);  % dxdxi
    
    x(ind) = X(m)+(xi+1).*Jacobian(m,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Righthandside

% f = exp(x).*(sin(pi*x)+pi*cos(pi*x));
c = 5;
f = 1-(1+c*x).*exp(c*x)/exp(c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology

D = spdiags([-ones(N*H,1) ones(N*H,1)],[0 1],N*H,N*H+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basis-functions

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass-matrix of zero-forms (flux-forms)

M0 = zeros(1,N*H+1);
for m=1:H
    ind = N*(m-1)+(1:N+1);
    M0(ind) = M0(ind) + w.*Jacobian(m,:);
end
M0 = diag(M0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass-matrix of one-forms (volume-forms)

M1 = zeros(H*N);
for m=1:H

    M1e = zeros(N);
    for i=1:N
        for j=1:N
            M1e(i,j) = sum(w.*e(i,:).*e(j,:)./Jacobian(m,:));
        end
    end
    ind = (m-1)*N+(1:N);
    M1(ind,ind) = M1e;
end
M1 = sparse(M1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix assembly and solve

L = [ M0      D'*M1
      M1*D  zeros(N*H) ];

RHS = [zeros(N*H+1,1) ; M1*diff(f)'];

% solve
qa = L\RHS;

q = qa(1:N*H+1);
a = qa(N*H+2:2*N*H+1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-proces

% Interpolation grid per element
kk = 100;
xixi=linspace(-1,1,kk);

[hh,ee]=MimeticpolyVal(xixi,N,1);

% Complete interpolation grid
XiXi = zeros(1,kk*H);
xx   = zeros(1,kk*H);
dxxdXiXi = zeros(1,kk*H);
for m=1:H
    ind = kk*(m-1)+(1:kk);
    XiXi(ind) = xixi+2*(m-1);
    xx(ind) = X(m)+(xixi+1)*(X(m+1)-X(m))/2;
    dxxdXiXi(ind(1:end-1)) = diff(xx(ind))./diff(XiXi(ind));
    dxxdXiXi(ind(end)) = dxxdXiXi(ind(end-1));
end

% Exact solution
% aa_ex = exp(xx).*sin(pi*xx);
aa_ex = xx.*(1-exp(c*xx)/exp(c));


% Numerical solution interpolated
aa = zeros(size(xx));
for m=1:H
    ind1 = kk*(m-1)+(1:kk);
    ind2 = N*(m-1)+(1:N);
    aa(ind1) = (a(ind2)'./Jacobian(m,1:N))*ee;
end

% Plot
plot(xx,aa_ex,'r')
hold on
plot(xx,aa,'--b')

for m=1:H
    for n=1:N
        i = n+(m-1)*N;
        plot([x(i) x(i+1)],[a(i) a(i)]/(x(n+1)-x(n)),'g')
    end
end
grid
