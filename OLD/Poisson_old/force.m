function [F] = force(N,x,y,c)

n = 8;

[xibar,Gw] = Gnodes(n); GGw = Gw'*Gw;
etabar = xibar;
F = zeros(N);
for i=1:N
    for j=1:N
        xi  = ( (x(i)+x(i+1))/2+(x(i+1)-x(i))/2*xibar )'*ones(1,n);
        eta = ones(n,1)*( (y(j)+y(j+1))/2+(y(j+1)-y(j))/2*etabar );

        xx = xi +c*sin(pi*xi).*sin(pi*eta);
        yy = eta+c*sin(pi*xi).*sin(pi*eta);
        
        f = -2*pi^2*sin(pi*xx).*sin(pi*yy);

        dxdxi  = 1+c*pi*cos(pi*xi).*sin(pi*eta);
        dxdeta = c*pi*sin(pi*xi).*cos(pi*eta);
        dydxi  = c*pi*cos(pi*xi).*sin(pi*eta);
        dydeta = 1+c*pi*sin(pi*xi).*cos(pi*eta);

        J = dxdxi.*dydeta-dxdeta.*dydxi;
        
        Jbar = ((x(i+1)-x(i))/2)*((y(j+1)-y(j))/2)*ones(n);

        F(i,j) = sum(sum(GGw.*f.*J.*Jbar));
        
    end
end

% sum(sum(F))