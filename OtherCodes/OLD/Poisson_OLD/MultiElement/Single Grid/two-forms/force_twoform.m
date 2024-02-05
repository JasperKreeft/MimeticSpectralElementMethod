function [F] = force_twoform(r,c)

global N numRows numColumns
global xi m cc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta = xi;

xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numRows+1)   ;   % This might be more advanced !!!!!

xiLR  = (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xi';
etaAB = (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*eta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 8;

[xin,Gw] = Gnodes(n); GGw = Gw'*Gw;
etan = xin;
F = zeros(N);
for i=1:N
    for j=1:N
        xibar  = ( (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xin )'*ones(1,n);
        etabar = ones(n,1)*( (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etan );

        xx = xibar +cc*sin(pi*xibar).*sin(pi*etabar);
        yy = etabar+cc*sin(pi*xibar).*sin(pi*etabar);

        f = -2*m^2*pi^2*sin(m*pi*xx).*sin(m*pi*yy);

        dxdxib  = 1+cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dxdetab = cc*pi*sin(pi*xibar).*cos(pi*etabar);
        dydxib  = cc*pi*cos(pi*xibar).*sin(pi*etabar);
        dydetab = 1+cc*pi*sin(pi*xibar).*cos(pi*etabar);

        dxibdxi   = (xiLR(i+1,1)-xiLR(i,1))/2*ones(n);
        dxibdeta  = zeros(n);
        detabdxi  = zeros(n);
        detabdeta = (etaAB(1,j+1)-etaAB(1,j))/2*ones(n);
        
        dxdxi  = dxdxib.*dxibdxi +dxdetab.*detabdxi ;
        dxdeta = dxdxib.*dxibdeta+dxdetab.*detabdeta;
        dydxi  = dydxib.*dxibdxi +dydetab.*detabdxi ;
        dydeta = dydxib.*dxibdeta+dydetab.*detabdeta;

        J = dxdxi.*dydeta-dxdeta.*dydxi;

%         J1 = (xiLR(i+1)-xiLR(i))/2 * (etaAB(j+1)-etaAB(j))/2;
%         
%         J2 = 1/numColumns * 1/numRows;

        F(i,j) = sum(sum(GGw.*f.*J));
    end
end