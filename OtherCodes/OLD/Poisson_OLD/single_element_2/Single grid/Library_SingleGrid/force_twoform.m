function [F] = force_twoform(x,y,FunctionType,Domain,DomInfo)

global N

n = 8;

[xibar,Gw] = Gnodes(n); GGw = Gw'*Gw;
etabar = xibar;
F = zeros(N);
for i=1:N
    for j=1:N
        xi  = ( (x(i)+x(i+1))/2+(x(i+1)-x(i))/2*xibar )'*ones(1,n);
        eta = ones(n,1)*( (y(j)+y(j+1))/2+(y(j+1)-y(j))/2*etabar );

        [xx,yy,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(xi,eta,Domain,DomInfo);

        f = exact_solution(xx,yy,FunctionType,'force');

        J = dxdxi.*dydeta-dxdeta.*dydxi;

        Jbar = ((x(i+1)-x(i))/2)*((y(j+1)-y(j))/2)*ones(n);

        F(i,j) = sum(sum(GGw.*f.*J.*Jbar));

    end
end
