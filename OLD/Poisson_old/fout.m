phi = zeros(N+2); l=0;
for i=1:N+2
    for j=1:N+2
        if i==1 || i==N+2
            phi(i,j) = 0;%phi_exact(i,j);
        elseif j==1 || j==N+2
            phi(i,j) = 0;%phi_exact(i,j);
        else            
        l=l+1;
        phi(i,j) = phi_in(l);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,Gw] = Gnodes(200); eta = xi;
Xi = xi'*ones(1,200); Eta = Xi';


X = Xi+c*sin(pi*Xi).*sin(pi*Eta);
Y = Eta+c*sin(pi*Xi).*sin(pi*Eta);

phi_ex = sin(m*pi*X).*sin(m*pi*Y);

phi_exact = sin(m*pi*XEGEG).*sin(m*pi*YEGEG);



hi = LagrangeVal(xi,N,3); hj = hi;
P = zeros(200);
for k=1:200
    for l=1:200
        for i=1:N
            for j=1:N
                P(k,l)=P(k,l)+phi(i+1,j+1)*hi(i+1,k)*hj(j+1,l);
            end
        end
    end
end


P_interp = zeros(200);
for k=1:200
    for l=1:200
        for i=1:N
            for j=1:N
                P_interp(k,l)=P_interp(k,l)+phi_exact(i+1,j+1)*hi(i+1,k)*hj(j+1,l);
            end
        end
    end
end

Wij = Gw'*Gw;

errorL2(N)        = sqrt( sum(sum( (P-phi_ex).^2.*Wij )) );
errorL2_interp(N) = sqrt( sum(sum( (P_interp-phi_ex).^2.*Wij )) );
