function PHI = potentialValue(FunctionType,Domain,DomInfo)

global N numElements numRows numColumns
global xi
global n xibar Gw

eta = xi;

% Main Mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = kron(ones(1,numRows+1),linspace(-1,1,numColumns+1));
etabAB = kron(linspace(-1,1,numRows+1),ones(1,numColumns+1));

GGw = Gw*Gw';
etabar = xibar;

PHI = zeros(N*N,numElements);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        
ind = c+(r-1)*(numColumns+1);
xiLR  = (xibLR(ind)+xibLR(ind+1))/2+(xibLR(ind+1)-xibLR(ind))/2*xi';
etaAB = (etabAB(ind)+etabAB(ind+numColumns+1))/2+(etabAB(ind+numColumns+1)-etabAB(ind))/2*eta;
        

for i=1:N
    for j=1:N
        xib  = ( (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xibar )*ones(1,n);
        etab = ones(n,1)*( (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etabar' );
        
        [xx,yy,dxdxib,dxdetab,dydxib,dydetab] = coordinatemap(xib,etab,Domain,DomInfo);

        phi = exact_solution(xx,yy,FunctionType,'two');

        dxibdxi   = (xiLR(i+1,1)-xiLR(i,1))/2*ones(n);
        dxibdeta  = zeros(n);
        detabdxi  = zeros(n);
        detabdeta = (etaAB(1,j+1)-etaAB(1,j))/2*ones(n);
        
        dxdxi  = dxdxib.*dxibdxi +dxdetab.*detabdxi ;
        dxdeta = dxdxib.*dxibdeta+dxdetab.*detabdeta;
        dydxi  = dydxib.*dxibdxi +dydetab.*detabdxi ;
        dydeta = dydxib.*dxibdeta+dydetab.*detabdeta;

        J = dxdxi.*dydeta-dxdeta.*dydxi;

        ij = i+(j-1)*N;
        PHI(ij,rc) = sum(sum(GGw.*phi.*J));
    end
end

    end
end


































%%%%%%%%%%%%%%%%%%%%%%%%%

% global N
% global n xibar Gw
% 
% etabar = xibar'; GGw = Gw*Gw';
% PHI_exact = zeros(N*N,1);
% for j=1:N
%     for i=1:N
%         xi  = ( (x(i)+x(i+1))/2+(x(i+1)-x(i))/2*xibar )*ones(1,n);
%         eta = ones(n,1)*( (y(j)+y(j+1))/2+(y(j+1)-y(j))/2*etabar );
% 
%         [xx,yy,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(xi,eta,Domain,DomInfo);
% 
%         p = exact_solution(xx,yy,FunctionType,'two');
% 
%         J = dxdxi.*dydeta-dxdeta.*dydxi;
% 
%         Jbar = ((x(i+1)-x(i))/2)*((y(j+1)-y(j))/2)*ones(n);
% 
%         ij = i+(j-1)*N;
%         PHI_exact(ij,1) = sum(sum(GGw.*p.*J.*Jbar));
%         
%     end
% end
