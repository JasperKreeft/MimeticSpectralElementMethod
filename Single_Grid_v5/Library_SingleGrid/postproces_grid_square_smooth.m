function [Meshp,hp,ep] = postproces_grid_square_smooth(Domain,DomInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights

global N numElements numRows numColumns
global nn

if isempty(numRows)
    numRows = sqrt(numElements);
    numColumns = numRows;
end

nn = 2*N; %N;%2*N;%
Meshp.nn = nn;
[Meshp.xip,Meshp.wp] = GLLnodes(nn-1); xip = Meshp.xip; etap = xip; wp = Meshp.wp;
% xip = linspace(-1,1,nn); etap = xip; [~,wp] = GLLnodes(nn-1); disp('!!! postprocess grid changed !!! Do NOT calculate errors')

Meshp.W = kron(wp,wp)';

[hp,ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Domain,'SinDeformGrid_01'); L=1; end
if strcmp(Domain,'SinDeformGrid'); L=2; end
if strcmp(Domain,'Vorticity'); L=DomInfo; end

z = [];
for h=1:numColumns
    z = [ z  (2*h-1)*L/(2*numColumns)-L/2+L/(2*numColumns)*xip' ];
end
Xi = repmat(z,nn,numRows);
z = [];
for h=1:numRows
    z = [ z  (2*h-1)*L/(2*numRows)-L/2+L/(2*numRows)*etap' ];
end
Eta = kron(z,ones(nn,numColumns));

[Meshp.X,Meshp.dXdXi]  = map(Xi,L,numColumns);
[Meshp.Y,Meshp.dYdEta] = map(Eta,L,numRows);

Meshp.dXdEta = zeros(nn^2,numElements);
Meshp.dYdXi  = zeros(nn^2,numElements);

Meshp.J = (Meshp.dXdXi.*Meshp.dYdEta-Meshp.dXdEta.*Meshp.dYdXi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%