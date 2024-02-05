function [A]=massmatrix(Ne,Ng,Pe,w,c,h,Le_j,Lo_j,xi,xj)


A = zeros(Ng);

for l=1:Ne

    p = Pe(l); % order of polynomial expansion

    k = p+1; % nr of integration points / order of the lobatto quadrature

    % Lobatto integration points
    z = xj(k+1,1:k+1);

    for i=1:p+1
        for j=1:p+1
            for m = 1:k+1
                if m==1
                    Pij(m) = 4*(i==1)*(j==1) * Lo_j(k+1,m,p)^2;
                elseif m==k+1
                    Pij(m) = 4*(i==p+1)*(j==p+1) * Lo_j(k+1,m,p)^2;
                else
                    Pij(m) = (1-z(m))^2 * (1+z(m))^2 / ((z(m)-xi(l,i)) * (z(m)-xi(l,j))) * Lo_j(k+1,m,p)^2;
                end
            end
            Theta(i,j) = 1/( p^2 * (p+1)^2 * Le_j(p+1,j) * Le_j(p+1,i) ) * Pij*w(k,1:k+1)';
        end
    end
    
    s = c(l,1):c(l,p+1);
    A(s,s) = A(s,s)+h(l)/2*Theta;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%