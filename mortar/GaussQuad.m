function [x,w] = GaussQuad(p,varargin)
%
%   GaussQuad determines the nodes and weights for gauss quadrature
%   
%   [x,w] = GaussQuad(p,varargin)
%
%   input:
%       p           :: number of cells
%       varargin    :: variable argument
%                      1) METHOD: allows the user to select which method to use
%                           -'GW'   : Golub-Welsch eigenvalue method (default for N<128)
%                                     best suited for when N is small. 
%                           -'FAST' : Glaser-Liu-Rokhlin fast algorithm
%                                     faster for large N    
%                      2) [A,B] : scales the nodes and weights for the finite interval [A,B]
%
%   output:
%       x           :: gauss nodes
%       w           :: gauss weights
%
%   Copyright 2004-2009 Chebfun
%   $Revision: 1.0 $  $Date: 10/05/2009 $
%   Revised by: Peter Kuystermans
%   $Revision: 2.0 $  $Date: 01/10/2011 $  

%-------------------------------------------------------------------------%
% input check                                                             %
%-------------------------------------------------------------------------%

    if p < 0
        error('Input should be a positive number');
    end

%-------------------------------------------------------------------------%
% compute number of points                                                %
%-------------------------------------------------------------------------%

    n = p+1;

%-------------------------------------------------------------------------%
% set defaults                                                            %
%-------------------------------------------------------------------------%

    interval = [-1,1]; 
    method = 'default'; 

%-------------------------------------------------------------------------%
% second input argument                                                   %
%-------------------------------------------------------------------------%

    if nargin > 1
        if isa(varargin{1},'double') && length(varargin{1}) == 2, interval = varargin{1}; end
        if isa(varargin{1},'char'), method = varargin{1}; end
        if length(varargin) == 2,
            if isa(varargin{2},'double') && length(varargin{2}) == 2, interval = varargin{2}; end
            if isa(varargin{2},'char'), method = varargin{2}; end
        end
    end

%-------------------------------------------------------------------------%
% decide on GW or FAST                                                    %
%-------------------------------------------------------------------------%

    %%% GW
    if (n < 128 || strcmpi(method,'GW')) && ~strcmpi(method,'fast') 
       m = n-1;
       beta = .5./sqrt(1-(2*(1:m)).^(-2));	% 3-term recurrence coefficientss
       T = diag(beta,1) + diag(beta,-1);  	% jacobi matrix
       [V,D] = eig(T);                     	% eigenvalue decomposition
       x = diag(D); [x,i] = sort(x);       	% legendre points
       w = 2*V(1,i).^2;                   	% weights
    %%% FAST   
    else                                                       
       [x ders] = alg0_Leg(n);            	% nodes and P_n'(x)
       w = 2./((1-x.^2).*ders.^2)';        	% weights
    end
    w = (2/sum(w))*w;                      	% normalise so that sum(w) = 2

%-------------------------------------------------------------------------%
% rescale interval                                                        %
%-------------------------------------------------------------------------%

    a = interval(1); b = interval(2);
    x = .5*( (x+1)*b - (x-1)*a);
    w = .5*(b-a)*w';

end

%-------------------------------------------------------------------------%
% FAST algorithm (functions)                                              %
%-------------------------------------------------------------------------%

%%% function: alg0_Leg   
%%% driver for FAST algorithm
function [roots ders] = alg0_Leg(n)     

    %%% compute coefficients of P_m(0), Pm'(0), m = 0,..,N
    Pm2 = 0; Pm1 = 1; Ppm2 = 0; Ppm1 = 0;
    for k = 0:n-1
        P = -k*Pm2/(k+1);
        Pp = ((2*k+1)*Pm1-k*Ppm2)/(k+1);
        Pm2 = Pm1; Pm1 = P; Ppm2 = Ppm1; Ppm1 = Pp;
    end

    roots = zeros(n,1); ders = zeros(n,1);          % allocate storage
    if mod(n,2),roots((n-1)/2) = 0; 
        ders((n+1)/2) = Pp;                         % zero is a root
    else
        [roots(n/2+1) ders(n/2+1)] = alg2_Leg(P,n); 
    end                                             % find first root
    [roots ders] = alg1_Leg(roots,ders);            % compute roots and derivatives
end

%-------------------------------------------------------------------------%

%%% function: alg1_Leg
%%% main algorithm
function [roots ders] = alg1_Leg(roots,ders) 

    n = length(roots);
    if mod(n,2), N = (n-1)/2; 
        s = 1;
    else
        N = n/2; s = 0; 
    end

    %%% approximate roots via asymptotic formula
    k = (n-2+s)/2:-1:1;
    theta = pi*(4*k-1)/(4*n+2);
    roots(((n+4-s)/2):end) = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

    m = 30;                                 % number of terms in taylor expansion
    u = zeros(1,m+1); up = zeros(1,m+1);
    for j = N+1:n-1
        x = roots(j);                       % previous root   
        h = roots(j+1) - x;             	% initial approx (via asymptotic formula)    
        M = 1/h;                           	% scaling

        %%% recurrence relation for legendre polynomials (scaled)
        u(1) = 0;   u(2) = ders(j)/M;  up(1) = u(2); up(m+1) = 0;
        for k = 0:m-2
            u(k+3) = (2*x*(k+1)/M*u(k+2)+(k-n*(n+1)/(k+1))*u(k+1)/M^2)./((1-x.^2)*(k+2));
            up(k+2) = (k+2)*u(k+3)*M;
        end

        %%% flip for more accuracy in inner product calculation
        u = u(m+1:-1:1);
        up = up(m+1:-1:1);

        %%% newton iteration
        hh = [ones(m,1) ; M];
        step = inf;
        l = 0;
        while (abs(step) > eps) && (l < 10)
            l = l + 1;
            step = (u*hh)/(up*hh);
            h = h - step;
            hh = [M;cumprod(M*h+zeros(m,1))];	% powers of h (fastest)
            hh = hh(end:-1:1);
        end

        %%% update
        roots(j+1) = x + h;
        ders(j+1) = up*hh;    
    end

    %%% nodes are symmetric
    roots(1:N+s) = -roots(n:-1:N+1);
    ders(1:N+s) = ders(n:-1:N+1);
    
end

%-------------------------------------------------------------------------%

%%% function: alg2_Leg
%%% find the first root (NOTE: P_n'(0)==0)
function [x1 d1] = alg2_Leg(Pn0,n) 

    %%% approximate first root via asymptotic formula
    k = ceil(n/2);
    theta = pi*(4*k-1)/(4*n+2);
    x1 = (1-(n-1)/(8*n^3)-1/(384*n^4)*(39-28./sin(theta).^2)).*cos(theta);

    m = 30;     % number of terms in Taylor expansion
    M = 1/x1; 	% scaling

    %%% recurrence relation for legendre polynomials
    u = zeros(1,m+1); up = zeros(1,m+1);
    u(1) = Pn0;
    for k = 0:2:m-2
        u(k+3) = (k-n*(n+1)/(k+1))*u(k+1)/(M^2*(k+2)); 
        up(k+2) = (k+2)*u(k+3)*M;
    end

    %%% flip for more accuracy in inner product calculation
    u = u(m+1:-1:1);
    up = up(m+1:-1:1);

    x1k = ones(m+1,1);

    %%% Newton iteration
    step = inf;
    l = 0;
    while (abs(step) > eps) && (l < 10)
        l = l + 1;
        step = (u*x1k)/(up*x1k);
        x1 = x1 - step;
        x1k = [1;cumprod(M*x1+zeros(m,1))];	% powers of h (fastest)
        x1k = x1k(end:-1:1);
    end
    d1 = up*x1k;

end

