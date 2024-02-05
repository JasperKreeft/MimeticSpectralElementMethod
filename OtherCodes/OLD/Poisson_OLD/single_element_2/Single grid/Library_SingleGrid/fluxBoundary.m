function Qbc = fluxBoundary(FunctionType,Domain,DomInfo,boundary)

global N xi

eta = xi;

n = 8;

[xibar,Gw] = Gnodes(n);
etabar = xibar;

arg1 = strcmp(boundary,'left');
arg2 = strcmp(boundary,'right');
arg3 = strcmp(boundary,'below');
arg4 = strcmp(boundary,'above');

if arg1 || arg2
    Qxi = zeros(1,N);
    for j = 1:N
        Xi  =  ((-1)*arg1+(+1)*arg2)*ones(1,n);
        Eta = (eta(j+1)+eta(j))/2+(eta(j+1)-eta(j))/2*etabar;

        [X,Y,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(Xi,Eta,Domain,DomInfo);

        [qx qy] = exact_solution(X,Y,FunctionType,'flux');

        qxi  = -dxdeta.*qy + dydeta.*qx;

        Jbar = (eta(j+1)-eta(j))/2;
        Qxi(1,j) = sum( qxi.*Gw*Jbar );
        Qbc = Qxi;
    end

elseif arg3 || arg4
    Qeta = zeros(N,1);
    for i = 1:N
        Xi  = (xi(i+1)+xi(i))/2+(xi(i+1)-xi(i))/2*xibar;
        Eta =  ((-1)*arg3+(+1)*arg4)*ones(1,n);

        [X,Y,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(Xi,Eta,Domain,DomInfo);

        [qx qy] = exact_solution(X,Y,FunctionType,'flux');

        qeta =   dxdxi.*qy -  dydxi.*qx;

        Jbar = ((xi(i+1)-xi(i))/2);
        Qeta(i,1) = sum( qeta.*Gw*Jbar);
        Qbc = Qeta;
    end
end