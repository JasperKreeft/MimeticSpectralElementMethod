function [x,w] = GaussQuad(p,varargin)
%LEGPTS  Legendre points and Gauss Quadrature Weights.
%  LEGPTS(N) returns N Legendre points in (-1,1).
%
%  [X,W] = LEGPTS(N) also returns a vector of weights for Gauss quadrature.
%
%  [X,W] = LEGPTS(N,METHOD) allows the user to select which method to use.
%       METHOD = 'GW' will use the traditional Golub-Welsch eigenvalue method,
%       which is best suited for when N is small. METHOD = 'FAST' will use 
%       the Glaser-Liu-Rokhlin fast algorithm, which is much faster for large N.
%       By default LEGPTS will use 'GW' when N < 128.
%
%  [X,W] = LEGPTS(N,[A,B]) scales the nodes and weights for the finite interval [A,B].
%
%  See http://www.comlab.ox.ac.uk/chebfun for chebfun information.

%  'GW' by Nick Trefethen, March 2009 - algoarithm adapted from [1].
%  'FAST' by Nick Hale, April 2009 - algorithm adapted from [2].
%
%  Copyright 2002-2009 by The Chebfun Team. 
%  Last commit: $Author: nich $: $Rev: 458 $:
%  $Date: 2009-05-10 20:51:03 +0100 (Sun, 10 May 2009) $:
%
%  References:
%   [1] G. H. Golub and J. A. Welsch, "Calculation of Gauss quadrature
%       rules", Math. Comp. 23:221-230, 1969, 
%   [2] A. Glaser, X. Liu and V. Rokhlin, "A fast algorithm for the 
%       calculation of the roots of special functions", SIAM Journal  
%       on Scientific Computing", 29(4):1420-1438:, 2007.

% compute the number of points
n = p+1;

if n <= 0
    error('Input should be a positive number');
end

% defaults
interval = [-1,1];
method = 'default';

% check inputs
if nargin > 1
    if isa(varargin{1},'double') && length(varargin{1}) == 2, interval = varargin{1}; end
    if isa(varargin{1},'char'), method = varargin{1}; end
    if length(varargin) == 2,
        if isa(varargin{2},'double') && length(varargin{2}) == 2, interval = varargin{2}; end
        if isa(varargin{2},'char'), method = varargin{2}; end
    end
end

% decide to use GW or FAST
if (n < 128 || strcmpi(method,'GW')) && ~strcmpi(method,'fast') % GW, see [1]
   m = n-1;
   beta = .5./sqrt(1-(2*(1:m)).^(-2));   % 3-term recurrence coeffs
   T = diag(beta,1) + diag(beta,-1);     % Jacobi matrix
   [V,D] = eig(T);                       % eigenvalue decomposition
   x = diag(D); [x,i] = sort(x);         % Legendre points
   w = 2*V(1,i).^2;                      % weights
else                                                            % Fast, see [2]
   [x ders] = alg0_Leg(n);               % nodes and P_n'(x)
   w = 2./((1-x.^2).*ders.^2)';          % weights
end
w = (2/sum(w))*w;                        % normalise so that sum(w) = 2

% rescale to arbitrary interval
a = interval(1); b = interval(2);
x = .5*( (x+1)*b - (x-1)*a);
w = .5*(b-a)*w';

% -------------------- Routines for FAST algorithm ------------------------

function [roots ders] = alg0_Leg(n) % driver for 'Fast'.
% Compute coefficients of P_m(0), Pm'(0), m = 0,..,N
Pm2 = 0; Pm1 = 1; Ppm2 = 0; Ppm1 = 0;
for k = 0:n-1
    P = -k*Pm2/(k+1);
    Pp = ((2*k+1)*Pm1-k*Ppm2)/(k+1);
    Pm2 = Pm1; Pm1 = P; Ppm2 = Ppm1; Ppm1 = Pp;
end

roots = zeros(n,1); ders = zeros(n,1);                      % allocate storage
if mod(n,2),roots((n-1)/2) = 0; ders((n+1)/2) = Pp;         % zero is a root
else [roots(n/2+1) ders(n/2+1)] = alg2_Leg(P,n); end        % find first root

[roots ders] = alg1_Leg(roots,ders); % compute roots and derivatives

% -------------------------------------------------------------------------

function [roots ders] = alg1_Leg(roots,ders) % main algorithm for 'Fast'
n = length(roots);
if mod(n,2), N = (n-1)/2; s = 1;
else N = n/2; s = 0; end

% Approximate roots via asymptotic formula
k = (n-2+s)/2:-1:1;
theta = pi*(4*k-1)/(4*n+2);
roots(((n+4-s)/2):end) = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

m = 30; % number of terms in Taylor expansion
u = zeros(1,m+1); up = zeros(1,m+1);
for j = N+1:n-1
    x = roots(j); % previous root
    
    % initial approx (via asymptotic foruma)
    h = roots(j+1) - x;
           
    % scaling
    M = 1/h;

    % recurrence relation for Legendre polynomials (scaled)
    u(1) = 0;   u(2) = ders(j)/M;  up(1) = u(2); up(m+1) = 0;
    for k = 0:m-2
        u(k+3) = (2*x*(k+1)/M*u(k+2)+(k-n*(n+1)/(k+1))*u(k+1)/M^2)./((1-x.^2)*(k+2));
        up(k+2) = (k+2)*u(k+3)*M;
    end
    
    % flip for more accuracy in inner product calculation
    u = u(m+1:-1:1);
    up = up(m+1:-1:1);
    
    hh = [ones(m,1) ; M];
    step = inf;
    l = 0;
    % Newton iteration
    while (abs(step) > eps) && (l < 10)
        l = l + 1;
        step = (u*hh)/(up*hh);
        h = h - step;
        hh = [M;cumprod(M*h+zeros(m,1))]; % powers of h (This is the fastest way!)
        hh = hh(end:-1:1);
    end
    
    % update
    roots(j+1) = x + h;
    ders(j+1) = up*hh;    
end

% nodes are symmetric
roots(1:N+s) = -roots(n:-1:N+1);
ders(1:N+s) = ders(n:-1:N+1);


% -------------------------------------------------------------------------

function [x1 d1] = alg2_Leg(Pn0,n) % find the first root (note P_n'(0)==0)

% % % advance ODE via Runge-Kutta for initial approx (redundant)
% x1 = rk2_Leg(0,-pi/2,0,n)

% Approximate first root via asymptotic formula
k = ceil(n/2);
theta = pi*(4*k-1)/(4*n+2);
x1 = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

m = 30; % number of terms in Taylor expansion

% scaling
M = 1/x1;
    
% recurrence relation for Legendre polynomials
u = zeros(1,m+1); up = zeros(1,m+1);
u(1) = Pn0;
for k = 0:2:m-2
    u(k+3) = (k-n*(n+1)/(k+1))*u(k+1)/(M^2*(k+2)); 
    up(k+2) = (k+2)*u(k+3)*M;
end

% flip for more accuracy in inner product calculation
u = u(m+1:-1:1);
up = up(m+1:-1:1);

x1k = ones(m+1,1);

step = inf;
l = 0;
% Newton iteration
while (abs(step) > eps) && (l < 10)
    l = l + 1;
    step = (u*x1k)/(up*x1k);
    x1 = x1 - step;
    x1k = [1;cumprod(M*x1+zeros(m,1))]; % powers of h (This is the fastest way!)
    x1k = x1k(end:-1:1);
end
d1 = up*x1k;


% -------------------------------------------------------------------------

% function x = rk2_Leg(t,tn,x,n) % Runge-Kutta for Legendre Equation (redundant)
% m = 10; h = (tn-t)/m;
% for j = 1:m
%     f1 = (1-x.^2);
%     k1 = -h*sqrt(f1)./(sqrt(n*(n+1)) - 0.5*x.*sin(2*t));
%     t = t+h;  f2 = (1-(x+k1).^2);
%     k2 = -h*sqrt(f2)./(sqrt(n*(n+1)) - 0.5*(x+k1).*sin(2*t));   
%     x = x+.5*(k1+k2);
% end