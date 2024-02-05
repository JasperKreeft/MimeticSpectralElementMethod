function [Qxi Qeta] = fluxValue(FunctionType,Domain,DomInfo)

% WERKT ALLEEN VOOR 1X1 ELEMENT

global N numRows numColumns 
global xi
global globalnr_1v globalnr_1h

eta = xi;

n = 8;

[xibar,Gw] = Gnodes(n);
etabar = xibar;

Qxi = zeros(size(globalnr_1v));
for i = 1:N+1
    for j = 1:N
        Xi  =  xi(i)*ones(1,n);
        Eta = (eta(j+1)+eta(j))/2+(eta(j+1)-eta(j))/2*etabar;

        [X,Y,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(Xi,Eta,Domain,DomInfo);

        [qx qy] = exact_solution(X,Y,FunctionType,'one');

        qxi  = -dxdeta.*qy + dydeta.*qx;

        Jbar = (eta(j+1)-eta(j))/2;
        Qxi(i,j) = sum( qxi.*Gw*Jbar );
    end
end

Qeta = zeros(size(globalnr_1h));
for i = 1:N
    for j = 1:N+1
        Xi  = (xi(i+1)+xi(i))/2+(xi(i+1)-xi(i))/2*xibar;
        Eta =  eta(j)*ones(1,n);
        
        [X,Y,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(Xi,Eta,Domain,DomInfo);

        [qx qy] = exact_solution(X,Y,FunctionType,'one');

        qeta =   dxdxi.*qy -  dydxi.*qx;
        
        Jbar = ((xi(i+1)-xi(i))/2);
        Qeta(i,j) = sum( qeta.*Gw*Jbar);
    end
end
