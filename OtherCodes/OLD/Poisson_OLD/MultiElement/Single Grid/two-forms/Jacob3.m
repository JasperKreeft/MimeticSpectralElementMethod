function Jb = Jacob3(xibar,etabar,W)

global cc numColumns numRows

% dxdxi  = 1+cc*pi*cos(pi*xibar).*sin(pi*etabar);
% dxdeta = cc*pi*sin(pi*xibar).*cos(pi*etabar);
% dydxi  = cc*pi*cos(pi*xibar).*sin(pi*etabar);
% dydeta = 1+cc*pi*sin(pi*xibar).*cos(pi*etabar);

        dxdxib  = 1+cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dxdetab = cc*pi*sin(pi*xibar).*cos(pi*etabar);
        dydxib  = cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dydetab = 1+cc*pi*sin(pi*xibar).*cos(pi*etabar);
        
        dxibdxi   = 1/numColumns;
        dxibdeta  = 0;
        detabdxi  = 0;
        detabdeta = 1/numRows;

        dxdxi  = dxdxib.*dxibdxi +dxdetab.*detabdxi ;
        dxdeta = dxdxib.*dxibdeta+dxdetab.*detabdeta;
        dydxi  = dydxib.*dxibdxi +dydetab.*detabdxi ;
        dydeta = dydxib.*dxibdeta+dydetab.*detabdeta;

J  = dxdxi.*dydeta-dxdeta.*dydxi;

Jb = sum(sum(W.*J))/4;
