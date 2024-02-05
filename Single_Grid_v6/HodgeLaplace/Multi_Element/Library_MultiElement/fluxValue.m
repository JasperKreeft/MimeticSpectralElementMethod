function [Qxi Qeta] = fluxValue(FunctionType,Domain,DomInfo)

global N numRows numColumns 
global xi
global globalnr_1v globalnr_1h

eta = xi;

n = 3*N;

[xibar,Gw] = GLnodes(n); Gw = Gw';
etabar = xibar;

xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numRows+1)   ;   % This might be more advanced !!!!!

Qxi = zeros(size(globalnr_1v));
Qeta = zeros(size(globalnr_1h));
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        
xiLR  = (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xi';
etaAB = (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*eta;


for i = 1:N+1
    for j = 1:N
        ij = i+(j-1)*(N+1);

        xib  = xiLR(i)*ones(n,1);
        etab = (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etabar;

        [xx,yy,dxdxib,dxdetab,dydxib,dydetab] = coordinatemap(xib,etab',Domain,DomInfo);

        [qx qy] = exact_solution(xx,yy,FunctionType,'one');

        dxibdeta  = zeros(n,1);
        detabdeta = (etaAB(1,j+1)-etaAB(1,j))/2*ones(n,1);

        dxdeta = dxdxib.*dxibdeta+dxdetab.*detabdeta;
        dydeta = dydxib.*dxibdeta+dydetab.*detabdeta;

        qxi  = -dxdeta.*qy + dydeta.*qx;

        Qxi(ij,rc) = sum( qxi.*Gw);

    end
end


for i = 1:N
    for j = 1:N+1
        ij = j+(i-1)*(N+1);

        xib  = ( (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xibar );
        etab = etaAB(j)*ones(n,1);

        [xx,yy,dxdxib,dxdetab,dydxib,dydetab] = coordinatemap(xib',etab,Domain,DomInfo);

        [qx qy] = exact_solution(xx,yy,FunctionType,'one');

        dxibdxi   = (xiLR(i+1,1)-xiLR(i,1))/2*ones(n,1);
        detabdxi  = zeros(n,1);
        
        dxdxi  = dxdxib.*dxibdxi +dxdetab.*detabdxi ;
        dydxi  = dydxib.*dxibdxi +dydetab.*detabdxi ;
        
        qeta =   dxdxi.*qy -  dydxi.*qx;

        Qeta(ij,rc) = sum( qeta.*Gw);
        
    end
end

    end
end