function [La,dLa,Laweights] = LagrangePoly(xi,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lagrange polynomials for order N-1 on a domain [min(xi) max(xi)],
% corresponding to a vector xi of length N.
%
% function [La,dLa] = LagrangePoly(xi,k)
%
% function returns a matrix La (size NxN) 
%
% if k==1 || isempty(k)==true
% then the rows of matrix La are the coefficients of the N-1 polynomials
% 
% if k==2
% then the first element equals 
% 1/((xi_i-xi_1)...(xi_i-xi_i-1)(xi_i-xi_i+1)...(xi_i-xi_N))
% and the other items in the row are zeros of the i-th Lagrange polynomial
% La(i,2:end) = [xi_1 ... xi_i-1 xi_i+1 ... xi_N]
%
% example:
% La = LagrangePoly([-1 0 1],1) gives
% 
% La = [  0.5000   -0.5000         0
%        -1.0000         0    1.0000
%         0.5000    0.5000         0 ]
% 
% La = LagrangePoly([-1 0 1],2) gives
% 
% La = [  0.5000         0    1.0000
%        -1.0000   -1.0000    1.0000
%         0.5000   -1.0000         0 ]
% 
% Written by Jasper Kreeft - 2009
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin==1
    k=1;
end

n = length(xi);

La = zeros(n);

ci = zeros(1,n);

for i=1:n
    xj     = xi;
    xj(i)  = [];
    ci(i)  = prod(xi(i)-xj);
    La(i,:) = [1/ci(i) xj];
end

LaPoly = zeros(n);
dLaPoly = zeros(n,n-1);
for i=1:n
    LaPoly(i,:) = La(i,1)*poly(La(i,2:n));
    dLaPoly(i,:) = polyder_J(LaPoly(i,:));
end
dLa = dLaPoly;

if k~=2
    La = LaPoly;
end

if nargout==3
    H = zeros(n,n+1);
    Laweights = zeros(1,n);
    for i=1:n
        H(i,:) = polyint(LaPoly(i,:));
        Laweights(i) = diff(polyval_J(H(i,:),[-1 1]));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To plot the Lagrange polynomials, uncommand the lines below

% kk = 100;
% x = linspace(min(xi),max(xi),kk);
% f = zeros(n,kk);
% for i=1:n
%     xj    = xi;
%     xj(i) = [];
%     for k=1:kk
%         f(i,k) = 1/ci(i)*prod(x(k)-xj);
%     end
% end
% 
% figure
% plot(x,f)
% grid
% legend(num2str((1:n)'))