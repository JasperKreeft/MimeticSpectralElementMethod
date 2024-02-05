function Gr = GrevilleAbscissa(X,P)

Xi = [ X(1)*ones(1,P) X X(end)*ones(1,P) ];

% Greville abscissa

n = length(X)+P-1;

Gr = zeros(1,n);

for j=1:n

    Gr(j) = sum(Xi(j+1:j+P))/P;

end