function L = reductionBspline(F,X,P)
% L = reductionBspline(F,X,P)
%
% Written by Jasper Kreeft, 2011
   
Xi = [ X(1)*ones(1,P) X X(end)*ones(1,P) ];

% Greville abscissa
Gr = GrevilleAbscissa(X,P);

% dual functionals
L = zeros(1,length(X)+P-1);
for j=1:length(X)+P-1

    xi_j_jp = Xi(j+1:j+P);

    n = length(xi_j_jp)+1;

    phi = zeros(P+1,n);
    phi(1,:) = poly(xi_j_jp);
    for a = 1:P
        phi(a+1,:) = [0 (n-1:-1:1).*phi(a,1:n-1)];
    end

    phiGr = zeros(P+1,1);
    for a=1:P+1
        phiGr(a) = phi(a,:)*(Gr(j).^(size(phi,2)-1:-1:0))';
    end


    for d=0:P
        L(j) = L(j) + (-1)^d/factorial(P)*phiGr(P-d+1)*F(d+1,j);
    end

end