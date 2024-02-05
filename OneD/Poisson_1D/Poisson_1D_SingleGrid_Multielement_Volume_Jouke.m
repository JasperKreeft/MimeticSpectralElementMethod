%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Simple 1D Poisson problem                                               %
% solved with Mimetic Spectral Element Method                             %
% Case: dd*a=f, where a is a 1-form (volume form)                         %
%                                                                         %
% Written by Jasper Kreeft - 2010                                         %
% Contact: j.j.kreeft@tudelft.nl                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings

H = 3; % number of spectral elements
N = 3; % number of cells in an element

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation

[xi,w] = GLLnodes(N);

Xi = zeros(1,N*H+1);
for m=1:H
    ind = N*(m-1)+(1:N+1);
    Xi(ind) = xi+2*(m-1);
end

Jacobian = 1/(2*H);

% Grid
x = (Xi+1)*Jacobian;

dxdXi = Jacobian;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Righthandside

% f = exp(x).*(sin(pi*x)+pi*cos(pi*x));
c = 15;
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
    M0(ind) = M0(ind) + w*dxdXi;
end
M0 = diag(M0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass-matrix of one-forms (volume-forms)

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(w.*e(i,:).*e(j,:)/dxdXi);
    end
end

M1 = kron(eye(H),M1); % Only possible for uniform grid !!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix assembly and solve

L = -M1*D/M0*D'*M1;

RHS = M1*diff(f)';
    
% solve
a = L\RHS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-proces

% Interpolation grid per element
kk = 100;
xixi=linspace(-1,1,kk);

[hh,ee]=MimeticpolyVal(xixi,N,1);

% Complete interpolation grid
XiXi = zeros(1,kk*H);
for m=1:H
    ind = kk*(m-1)+(1:kk);
    XiXi(ind) = xixi+2*(m-1);
end
xx = (XiXi+1)*Jacobian;

% Exact solution
% aa_ex = exp(xx).*sin(pi*xx);
aa_ex = xx.*(1-exp(c*xx)/exp(c));

% Numerical solution interpolated
aa = zeros(size(xx));
for m=1:H
    ind1 = kk*(m-1)+(1:kk);
    ind2 = N*(m-1)+(1:N);
    aa(ind1) = a(ind2)'*ee/Jacobian;
end

% Plot
plot(xx,aa_ex,'r')
hold on
plot(xx,aa,'--b')
for i=1:N*H
    plot([x(i) x(i+1)],[a(i) a(i)]/(x(i+1)-x(i)),'g')
end
grid
