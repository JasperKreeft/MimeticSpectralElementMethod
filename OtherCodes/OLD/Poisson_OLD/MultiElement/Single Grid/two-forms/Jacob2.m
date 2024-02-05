function JJ = Jacob2(xibar,etabar)

global cc


dxdxi  = 1+cc*pi*cos(pi*xibar).*sin(pi*etabar);
dxdeta = cc*pi*sin(pi*xibar).*cos(pi*etabar);
dydxi  = cc*pi*cos(pi*xibar).*sin(pi*etabar);
dydeta = 1+cc*pi*sin(pi*xibar).*cos(pi*etabar);

JJ  = dxdxi.*dydeta-dxdeta.*dydxi;
        
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


