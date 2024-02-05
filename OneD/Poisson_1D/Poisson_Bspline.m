clear all
close all
clc

path(path,'/media/My Passport/matlabfunctions/IsoRene')

X = [-1 1];
Nx = length(X);

kk = 1000;
[xx,ww] = GLLnodes(kk);

Prange = 2:15;

L2 = zeros(size(Prange));
C  = zeros(size(Prange));
n=0;
for P=Prange
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
n=n+1;
% L2(n) = sqrt( 2/kk*sum((aa-aa_ex).^2) );
L2(n) = sqrt( sum( (aa-aa_ex).^2.*ww ) );
C(n) = cond(Matrix);
% pause

end

nhf = floor(n/2);

figure
semilogy(Prange,L2,'-')
title('L_2 convergence')
grid on

figure
Line2 = 0.05*Prange.^3;
semilogy(Prange,Line2,'--g')
str2 = 'y=0.05x^3';
hold on
rate = (log(C(n))-log(C(nhf)))/(log(Prange(n))-log(Prange(nhf)));
% Line = C(n)/Prange(n)^rate*Prange.^rate;
% semilogy(Prange,Line,'r')
str = ['y=' num2str(C(n)/Prange(n)^rate,3) 'x^{' num2str(rate,3) '}'];
semilogy(Prange,C,'-');
grid on
y=ylim;
xlim([Prange(1) Prange(n)])
ylim([1 y(2)])
legend(str2,'Condition number',4)
% legend(str2,'Condition number',str,4)