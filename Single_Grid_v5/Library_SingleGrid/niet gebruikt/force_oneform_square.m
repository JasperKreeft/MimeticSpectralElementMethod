function [Fxi Feta] = force_oneform_square(i,FunctionType,Domain,DomInfo)

global N numElements numRows numColumns
global xi
global n xibar Gw

if isempty(numElements)
    numElements = 1;
end
if ~exist('numRows','var')
    numRows    = sqrt(numElements);
    numColumns = numRows;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta = xi;

xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numRows+1)   ;   % This might be more advanced !!!!!

c = i-floor((i-1)/numColumns)*numColumns;
r = ceil(i/numColumns);

xiLR  = (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xi';
etaAB = (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*eta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

etabar = xibar;

Fxi = zeros(N*(N+1),1);
for j = 1:N
    for i = 1:N+1
        ij = i+(j-1)*(N+1);

        xib  = xiLR(i)*ones(n,1);
        etab = (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etabar;

        [xx,yy,dxdxib,dxdetab,dydxib,dydetab] = coordinatemap(xib,etab,Domain,DomInfo);

        [fx,fy] = exact_solution(xx,yy,FunctionType,'force');

        dxibdeta  = zeros(n,1);
        detabdeta = (etaAB(1,j+1)-etaAB(1,j))/2*ones(n,1);
        
        dxdeta = dxdxib.*dxibdeta+dxdetab.*detabdeta;
        dydeta = dydxib.*dxibdeta+dydetab.*detabdeta;

        fxi  = -dxdeta.*fy + dydeta.*fx;

        Fxi(ij,1) = sum( fxi.*Gw);
    end
end

Feta = zeros(N*(N+1),1);
for i = 1:N
    for j = 1:N+1
        ij = j+(i-1)*(N+1);

        xib  = ( (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xibar );
        etab = etaAB(j)*ones(n,1);

        [xx,yy,dxdxib,dxdetab,dydxib,dydetab] = coordinatemap(xib,etab,Domain,DomInfo);

        [fx,fy] = exact_solution(xx,yy,FunctionType,'force');

        dxibdxi   = (xiLR(i+1,1)-xiLR(i,1))/2*ones(n,1);
        detabdxi  = zeros(n,1);
        
        dxdxi  = dxdxib.*dxibdxi +dxdetab.*detabdxi ;
        dydxi  = dydxib.*dxibdxi +dydetab.*detabdxi ;
        
        feta =   dxdxi.*fy -  dydxi.*fx;

        Feta(ij,1) = sum( feta.*Gw);

    end
end
