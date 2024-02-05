function U = u_exact(xiLR,etaAB)

global N m cc

n = 8;

[xin,Gw] = Gnodes(n); GGw = Gw'*Gw;
etan = xin;

% Jb = (xiLR(N+1)-xiLR(1))/2 * (etaAB(N+1)-etaAB(1))/2;

U = zeros(N);
for i=1:N
    for j=1:N
        xibar  = ( (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xin )'*ones(1,n);
        etabar = ones(n,1)*( (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etan );

        xx = xibar +cc*sin(pi*xibar).*sin(pi*etabar);
        yy = etabar+cc*sin(pi*xibar).*sin(pi*etabar);

        u = sin(m*pi*xx).*sin(m*pi*yy);

        dxdxib  = 1+cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dxdetab = cc*pi*sin(pi*xibar).*cos(pi*etabar);
        dydxib  = cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dydetab = 1+cc*pi*sin(pi*xibar).*cos(pi*etabar);

        dxibdxi   = (xiLR(i+1)-xiLR(i))/2;
        dxibdeta  = 0;
        detabdxi  = 0;
        detabdeta = (etaAB(j+1)-etaAB(j))/2;

        dxdxi  = dxdxib.*dxibdxi +dxdetab.*detabdxi ;
        dxdeta = dxdxib.*dxibdeta+dxdetab.*detabdeta;
        dydxi  = dydxib.*dxibdxi +dydetab.*detabdxi ;
        dydeta = dydxib.*dxibdeta+dydetab.*detabdeta;

        J  = dxdxi.*dydeta-dxdeta.*dydxi;

        U(i,j) = sum(sum(GGw.*u.*J));

    end
end
