function PHI_exact = potentialValue(x,y,FunctionType,Domain,DomInfo)

global N
global n xibar Gw

etabar = xibar'; GGw = Gw*Gw';
PHI_exact = zeros(N*N,1);
for j=1:N
    for i=1:N
        xi  = ( (x(i)+x(i+1))/2+(x(i+1)-x(i))/2*xibar )*ones(1,n);
        eta = ones(n,1)*( (y(j)+y(j+1))/2+(y(j+1)-y(j))/2*etabar );

        [xx,yy,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(xi,eta,Domain,DomInfo);

        p = exact_solution(xx,yy,FunctionType,'two');

        J = dxdxi.*dydeta-dxdeta.*dydxi;

        Jbar = ((x(i+1)-x(i))/2)*((y(j+1)-y(j))/2)*ones(n);

        ij = i+(j-1)*N;
        PHI_exact(ij,1) = sum(sum(GGw.*p.*J.*Jbar));
        
    end
end
