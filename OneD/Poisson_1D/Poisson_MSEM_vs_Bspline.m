clear all
close all
clc

Prange = 2:15;

kk = 1000;
[xx,ww] = GLLnodes(kk);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Mimetic Spectral Elements following the ideas of                      %
%   Kreeft, Palha and Gerritsma                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L2  = zeros(size(Prange));
C   = zeros(size(Prange));
DoF = zeros(size(Prange));
n=0;
for P=Prange

n=n+1;

N = P;

[xgl,wgl] = GLLnodes(N);

Jac = 1;

[h,e]=MimeticpolyVal(xgl,N,1);

NG = topology1D(N);

M0 = diag(wgl*Jac);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(wgl.*e(i,:).*e(j,:))/Jac;
    end
end

f = cos(pi/2*xgl);

Matrix = -NG'*M1*NG;
Rhs = M0*f';

%% boundary conditions
Matrix_full = Matrix;

Matrix(N+1,:) = [];
Matrix(1,:)   = [];
Matrix(:,N+1) = [];
Matrix(:,1)   = [];

Rhs_full = Rhs;
Rhs = Rhs(2:N);

% solve
a = Matrix\Rhs;

DoF(n) = length(a);

a = [ 0  a' 0 ];

%%

[hh,ee]=MimeticpolyVal(xx,N,1);

aa    = a*hh;
aa_ex = -4/pi^2*cos(pi/2*xx);

if length(Prange)==1
figure
plot(xx,aa_ex,'r')
hold on
plot(xx,aa,'--b')
plot(xgl,a,'x')
grid
end

%%

L2(n) = sqrt( sum( (aa-aa_ex).^2.*ww ) );
C(n) = cond(Matrix);

end

nhf = floor(n/2);

figure(1)
subplot(1,2,1)
semilogy(Prange,L2,'-r')
title('L_2 convergence vs P')
grid on
hold on
subplot(1,2,2)
semilogy(DoF,L2,'-r')
title('L_2 convergence vs DoF')
grid on
hold on

figure(2)
subplot(1,2,1)
loglog(Prange,C,'-r');
grid on
hold on
xlim([Prange(1) Prange(n)])
subplot(1,2,2)
semilogy(Prange,C,'-r');
grid on
hold on
xlim([Prange(1) Prange(n)])




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   B-splines following the ideas of Buffa, Evans and Hiemstra            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


addpath('c:\Users\Jasper.Kreeft\Data\MATLAB\Toolboxes\Bsplines')

X = [-1 1];
Nx = length(X);

L2  = zeros(size(Prange));
C   = zeros(size(Prange));
DoF = zeros(size(Prange));
n=0;
for P=Prange

n=n+1;

N=P;

[xgl,wgl] = GLLnodes(N);

[B,E]=Bspline(xgl,X,P);
[Bb,Ee]=Bspline(xx,X,P);

NG = topology1D(P+Nx-2);

M0 = zeros(P+Nx-1);
for i=1:P+Nx-1
    for j=1:P+Nx-1
        M0(i,j) = sum(wgl.*B(i,:).*B(j,:));
    end
end

M1 = zeros(P+Nx-2);
for i=1:P+Nx-2
    for j=1:P+Nx-2
        M1(i,j) = sum(wgl.*E(i,:).*E(j,:));
    end
end

%% f = sin(pi*xgl);

% Greville abscissa
Gr = GrevilleAbscissa(X,P);

% Function and its derivatives

% F(1,:) = cos(pi/2*Gr);
F = zeros(P+1,length(Gr));
for d=0:P
F(d+1,:) = (-1)^ceil(d/2)*(pi/2)^d*(even(d)*cos(pi/2*Gr)+odd(d)*sin(pi/2*Gr));
end

% Reduction of function using dual functionals
f = reductionBspline(F,X,P);

if length(Prange)==1
fh = f*Bb;
subplot(1,2,1)
plot(xx,fh,'g')
end

%%

Matrix = -NG'*M1*NG;
Rhs = M0*f';

%% boundary conditions
Matrix_full = Matrix;

Matrix(P+1,:) = [];
Matrix(1,:)   = [];
Matrix(:,P+1) = [];
Matrix(:,1)   = [];

Rhs_full = Rhs;
Rhs = Rhs(2:P+Nx-2);

% solve
a = Matrix\Rhs;

DoF(n) = length(a);

a = [ 0  a' 0 ];

%%

aa    = a*Bb;
aa_ex = -4/pi^2*cos(pi/2*xx);

if length(Prange)==1
subplot(1,2,2)
plot(xx,aa_ex,'r')
hold on
plot(xx,aa,'--b')
plot(Gr,a,'x')
grid
end

%%

L2(n) = sqrt( sum( (aa-aa_ex).^2.*ww ) );
C(n) = cond(Matrix);

end

nhf = floor(n/2);

figure(1)
subplot(1,2,1)
semilogy(Prange,L2,'-')
title('L_2 convergence vs P')
grid on
hold on
subplot(1,2,2)
semilogy(DoF,L2,'-')
title('L_2 convergence vs DoF')
grid on
hold on


figure(2)
subplot(1,2,1)
loglog(Prange,C,'-');
grid on
y=ylim;
xlim([Prange(1) Prange(n)])
% ylim([1 y(2)])
subplot(1,2,2)
semilogy(Prange,C,'-');
grid on
y=ylim;
xlim([Prange(1) Prange(n)])
% ylim([1 y(2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(Prange)>1
figure(1)
subplot(1,2,1)
legend('MSEM','B-splines')
subplot(1,2,2)
legend('MSEM','B-splines')

figure(2)
subplot(1,2,1)
legend('MSEM','B-splines')
subplot(1,2,2)
legend('MSEM','B-splines')
end
