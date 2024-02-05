% hymanplot;
figure; axis([-1.5 1.5 -1.5 1.5]); hold on
clear all; clc;

N = 16;

xp = GLLnodes(N);
xd = Gnodes(N);
xd_ex = [-1 xd 1];
yp = xp; yd = xd; yd_ex = xd_ex;

c = 0.6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 50;

[xibar,Gw] = Gnodes(n); GGw = Gw'*Gw;
etabar = xibar;
F = zeros(N);
for i=1:N
    for j=1:N
        xi  = ( (xp(i)+xp(i+1))/2+(xp(i+1)-xp(i))/2*xibar )'*ones(1,n);
        eta = ones(n,1)*( (yp(j)+yp(j+1))/2+(yp(j+1)-yp(j))/2*etabar );

        xx = xi +c*sin(pi*xi).*sin(pi*eta);
        yy = eta+c*sin(pi*xi).*sin(pi*eta);

        f = sin(pi*xx).*sin(pi*yy); %-2*pi^2*

        figure(1)
%         subplot(1,2,2)
%         contourf(xx,yy,f)%surf(xx,yy,f)%
%         plot(xx,yy,'.k')
        pcolor(xx,yy,f); shading interp

        dxdxi  = 1+c*pi*cos(pi*xi).*sin(pi*eta);
        dxdeta = c*pi*sin(pi*xi).*cos(pi*eta);
        dydxi  = c*pi*cos(pi*xi).*sin(pi*eta);
        dydeta = 1+c*pi*sin(pi*xi).*cos(pi*eta);

        J = dxdxi.*dydeta-dxdeta.*dydxi;

        Jbar = ((xp(i+1)-xp(i))/2)*((yp(j+1)-yp(j))/2)*ones(n);

        F(i,j) = sum(sum(GGw.*f.*J.*Jbar));

    end
end