function [Ng,xe,xg,xen,xi,xj,h,c,Le,Le_j,w,Lo_j] = hp(Ne,Pe,L)

maxPe = max(Pe);
Ng = sum(Pe+1)-(Ne-1);    % number of global nodes

xe = linspace(0,1,Ne+1);
h  = diff(xe);

% Gauss-Lobatto zeros, element node and global node locations
% and connectivity matrix
P11 = JacobiPoly(maxPe,1,1); % p one higher than needed for xi. Needed for quadrature points z
xi = zeros(Ne,maxPe+1);
xen = zeros(Ne,maxPe+1);
xg = [];
c = zeros(Ne,maxPe+1);
Ic = 2;
for l=1:Ne
    p = Pe(l);
keyboard
    % Gauss-Lobatto-Legendre zeros
    innerzeros = roots(P11(p,p:-1:1))';
    xi(l,1:p+1) = sort([-1 innerzeros 1]);

    % Element nodes
    xen(l,1:p+1) = (xe(l)+xe(l+1))/2+xi(l,1:p+1)*(xe(l+1)-xe(l))/2;

    % Global nodes
    xg = [xg xen(l,1:p)];

    % Connectivity matrix
    Ic = Ic-1;
    for j=1:p+1
        c(l,j) = Ic;
        Ic = Ic+1;
    end

end
xg = [xg L];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomials                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Legendre polynomials
Le = JacobiPoly(maxPe+1,0,0);


% Gauss-Lobatto zeros
p = 0;
xj(1,1) = 1;
Le_j(1,1) = Le(1,1)*xj(1,1);
for p=1:maxPe+1
    innerzeros = roots(P11(p,p:-1:1))';
    xj(p+1,1:p+1) = sort([-1 innerzeros 1]);
    for j=1:p+1
        Le_j(p+1,j) = sum(Le(p+1,1:p+1).*xj(p+1,j).^(0:p));
    end
end

% Corresponding weights
w = zeros(maxPe+1,maxPe+2);
for p=0:maxPe
    % number of quadrature points
    k=p+1;
    for j=1:k+1
        w(k,j) = 2/(k*(k+1)*(Le_j(k+1,j))^2);
    end
end

% Lobatto polynomial
Lo_j = zeros(maxPe);
for p=1:maxPe+1
    Lo_pm1 = (p+1)/2*P11(p,1:p);
    for j=1:maxPe+2
        for k=1:j
            Lo_j(j,k,p) = sum(Lo_pm1.*xj(j,k).^(0:p-1));
        end
    end
end
