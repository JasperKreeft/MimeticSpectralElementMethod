function varargout = massfunction(varargin)

global N numElements numRows numColumns
global xi
global n xibar Gw

if isempty(numElements); numElements = 1; end
if ~exist('numRows','var'); numRows = sqrt(numElements); numColumns = numRows; end

eta = xi;
        
r            = varargin{1};
c            = varargin{2};
FunctionType = varargin{3};
Domain       = varargin{4};
DomInfo      = varargin{5};

% Main Mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = kron(ones(1,numRows+1),linspace(-1,1,numColumns+1));
etabAB = kron(linspace(-1,1,numRows+1),ones(1,numColumns+1));

ind = c+(r-1)*(numColumns+1);
xiLR  = (xibLR(ind)+xibLR(ind+1))/2+(xibLR(ind+1)-xibLR(ind))/2*xi';
etaAB = (etabAB(ind)+etabAB(ind+numColumns+1))/2+(etabAB(ind+numColumns+1)-etabAB(ind))/2*eta;

GGw = Gw*Gw';
etabar = xibar;

G = zeros(N*N,1);
for i=1:N
    for j=1:N
        xib  = ( (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xibar )*ones(1,n);
        etab = ones(n,1)*( (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etabar' );
        
        [xx,yy,dxdxib,dxdetab,dydxib,dydetab] = coordinatemap(xib,etab,Domain,DomInfo);

        g = exact_solution(xx,yy,FunctionType,'mass');

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
        G(ij,1) = sum(sum(GGw.*g.*J));
    end
end

varargout{1} = G;