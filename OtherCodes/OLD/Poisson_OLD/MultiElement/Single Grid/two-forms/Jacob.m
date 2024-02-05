function Jb = Jacob(xiLR,etaAB)

global numColumns numRows cc

n = 8;

[xin,Gw] = Gnodes(n); GGw = Gw'*Gw;
etan = xin;

Jb = zeros(numColumns,numRows);
for c=1:numColumns
    for r=1:numRows
        xibar  = ( (xiLR(c)+xiLR(c+1))/2+(xiLR(c+1)-xiLR(c))/2*xin )'*ones(1,n);
        etabar = ones(n,1)*( (etaAB(r)+etaAB(r+1))/2+(etaAB(r+1)-etaAB(r))/2*etan );

        dxdxi  = 1+cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dxdeta = cc*pi*sin(pi*xibar).*cos(pi*etabar);
        dydxi  = cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dydeta = 1+cc*pi*sin(pi*xibar).*cos(pi*etabar);

%         dxdxib  = 1+cc*pi*cos(pi*xibar).*sin(pi*etabar);
%         dxdetab = cc*pi*sin(pi*xibar).*cos(pi*etabar);
%         dydxib  = cc*pi*cos(pi*xibar).*sin(pi*etabar);
%         dydetab = 1+cc*pi*sin(pi*xibar).*cos(pi*etabar);
% 
%         dxibdxi   = (xiLR(c+1)-xiLR(c))/2;
%         dxibdeta  = 0;
%         detabdxi  = 0;
%         detabdeta = (etaAB(r+1)-etaAB(r))/2;
% 
%         dxdxi  = dxdxib.*dxibdxi +dxdetab.*detabdxi ;
%         dxdeta = dxdxib.*dxibdeta+dxdetab.*detabdeta;
%         dydxi  = dydxib.*dxibdxi +dydetab.*detabdxi ;
%         dydeta = dydxib.*dxibdeta+dydetab.*detabdeta;

        J  = dxdxi.*dydeta-dxdeta.*dydxi;

        Jb(c,r) = sum(sum(GGw.*J));

    end
end