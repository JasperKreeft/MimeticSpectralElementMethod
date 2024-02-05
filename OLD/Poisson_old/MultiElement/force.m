function [F] = force(x,y)

global N m cc

n = 8;

[xibar,Gw] = Gnodes(n); GGw = Gw'*Gw;
etabar = xibar;
F = zeros(N);
for i=1:N
    for j=1:N
        xi  = ( (x(i)+x(i+1))/2+(x(i+1)-x(i))/2*xibar )'*ones(1,n);
        eta = ones(n,1)*( (y(j)+y(j+1))/2+(y(j+1)-y(j))/2*etabar );

        xx = xi +cc*sin(pi*xi).*sin(pi*eta);
        yy = eta+cc*sin(pi*xi).*sin(pi*eta);

        f = -2*m^2*pi^2.*sin(m*pi*xx).*sin(m*pi*yy);

        dxdxi  = 1+cc*pi*cos(pi*xi).*sin(pi*eta);
        dxdeta = cc*pi*sin(pi*xi).*cos(pi*eta);
        dydxi  = cc*pi*cos(pi*xi).*sin(pi*eta);
        dydeta = 1+cc*pi*sin(pi*xi).*cos(pi*eta);

        J = dxdxi.*dydeta-dxdeta.*dydxi;
        
        Jbar = ((x(i+1)-x(i))/2)*((y(j+1)-y(j))/2)*ones(n);

        F(i,j) = sum(sum(GGw.*f.*J.*Jbar));
        
    end
end

% sum(sum(F))