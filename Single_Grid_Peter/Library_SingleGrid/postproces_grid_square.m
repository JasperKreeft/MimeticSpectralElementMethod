function [Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights

global N numElements numRows numColumns
global nn

if isempty(numRows)
    numRows = sqrt(numElements);
    numColumns = numRows;
end

nn = 4*N; % 
[xip,wp] = GLLnodes(nn-1); etap = xip;
Meshp.W = kron(wp,wp)';

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = kron(ones(1,numRows+1),linspace(-1,1,numColumns+1));
etabAB = kron(linspace(-1,1,numRows+1),ones(1,numColumns+1));

% xibLR  = kron(ones(1,numRows+1),GLLnodes(numColumns));
% etabAB = kron(GLLnodes(numRows),ones(1,numColumns+1));


dXibdXi   = zeros(nn^2,numElements);
dXibdEta  = zeros(nn^2,numElements);
dEtabdXi  = zeros(nn^2,numElements);
dEtabdEta = zeros(nn^2,numElements);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        ind = c+(r-1)*(numColumns+1);
        dXibdXi(:,rc)   = (xibLR(ind+1)-xibLR(ind))/2*ones(nn^2,1);
        dXibdEta(:,rc)  = zeros(nn^2,1);
        dEtabdXi(:,rc)  = zeros(nn^2,1);
        dEtabdEta(:,rc) = (etabAB(ind+numColumns+1)-etabAB(ind))/2*ones(nn^2,1);
    end
end

Etarc = zeros(nn^2,numElements);
Xirc = zeros(nn^2,numElements);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

        ind = c+(r-1)*(numColumns+1);
        Xirc(:,rc) = (xibLR(ind+1)+xibLR(ind))/2+(xibLR(ind+1)-xibLR(ind))/2*kron(ones(nn,1),xip');
        Etarc(:,rc) = (etabAB(ind+numColumns+1)+etabAB(ind))/2+(etabAB(ind+numColumns+1)-etabAB(ind))/2*kron(etap',ones(nn,1));

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