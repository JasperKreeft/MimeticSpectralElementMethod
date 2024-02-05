function [F] = force(xiLR,etaAB)

global N m cc

n = 8;

[xi,Gw] = Gnodes(n); GGw = Gw'*Gw;
eta = xi;
F = zeros(N);
for i=1:N
    for j=1:N
        xibar  = ( (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xi )'*ones(1,n);
        etabar = ones(n,1)*( (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*eta );

        xx = xibar +cc*sin(pi*xibar).*sin(pi*etabar);
        yy = etabar+cc*sin(pi*xibar).*sin(pi*etabar);

        f = -2*m^2*pi^2.*sin(m*pi*xx).*sin(m*pi*yy);

        dxdxib  = 1+cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dxdetab = cc*pi*sin(pi*xibar).*cos(pi*etabar);
        dydxib  = cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dydetab = 1+cc*pi*sin(pi*xibar).*cos(pi*etabar);

        dxibdxi   = (xiLR(i+1,1)-xiLR(i,1))/2*ones(n);
        dxibdeta  = zeros(n);
        detabdxi  = zeros(n);
        detabdeta = (etaAB(1,j+1)-etaAB(1,j))/2*ones(n);

        dxdxi  = dxdxib.*dxibdxi +dxdetab.*detabdxi ;
        dxdeta = dxdxib.*dxibdeta+dxdetab.*detabdeta;
        dydxi  = dydxib.*dxibdxi +dydetab.*detabdxi ;
        dydeta = dydxib.*dxibdeta+dydetab.*detabdeta;

        J = dxdxi.*dydeta-dxdeta.*dydxi;

        F(i,j) = sum(sum(GGw.*f.*J));
    end
end