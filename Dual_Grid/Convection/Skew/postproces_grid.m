function [Meshp,hp,ep] = postproces_grid(Domain,DomInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 

global N numColumns numRows
global nn

nn = 4*N; % 50;
[xip,wp] = GLLnodes(nn-1); etap = xip;
Meshp.W = wp'*wp;

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xibLR  = linspace(-1,1,numColumns+1);
etabAB = linspace(-1,1,numRows+1)   ;


dXibdXi   = zeros(nn,nn,numRows*numColumns);
dXibdEta  = zeros(nn,nn,numRows*numColumns);
dEtabdXi  = zeros(nn,nn,numRows*numColumns);
dEtabdEta = zeros(nn,nn,numRows*numColumns);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        dXibdXi(:,:,rc)   = (xibLR(c+1)-xibLR(c))/2*ones(nn);
        dXibdEta(:,:,rc)  = zeros(nn);
        dEtabdXi(:,:,rc)  = zeros(nn);
        dEtabdEta(:,:,rc) = (etabAB(r+1)-etabAB(r))/2*ones(nn);
    end
end

Etarc = zeros(nn,nn,numRows*numColumns);
Xirc = zeros(nn,nn,numRows*numColumns);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

        Xirc(:,:,rc) = ((xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xip)'*ones(1,nn);
        Etarc(:,:,rc) = ones(nn,1)*((etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*etap);
        
    end
end

[Xp,Yp,dXdXib,dXdEtab,dYdXib,dYdEtab] = coordinatemap(Xirc,Etarc,Domain,DomInfo);

Meshp.X = Xp;
Meshp.Y = Yp;

Meshp.dXdXi  = dXdXib.*dXibdXi +dXdEtab.*dEtabdXi ;
Meshp.dXdEta = dXdXib.*dXibdEta+dXdEtab.*dEtabdEta;
Meshp.dYdXi  = dYdXib.*dXibdXi +dYdEtab.*dEtabdXi ;
Meshp.dYdEta = dYdXib.*dXibdEta+dYdEtab.*dEtabdEta;

Meshp.J = (Meshp.dXdXi.*Meshp.dYdEta-Meshp.dXdEta.*Meshp.dYdXi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%