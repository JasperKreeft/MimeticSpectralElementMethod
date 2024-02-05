function varargout = forcefunction(form,varargin)


global N numElements numRows numColumns
global xi
global n xibar Gw

if isempty(numElements); numElements = 1; end
if ~exist('numRows','var'); numRows = sqrt(numElements); numColumns = numRows; end

eta = xi;


switch form

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 0
        
X            = varargin{1};
Y            = varargin{2};
FunctionType = varargin{3};
        
F = exact_solution(X,Y,FunctionType,'force');

varargout{1} = F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 1
        
r            = varargin{1};
c            = varargin{2};
FunctionType = varargin{3};
Domain       = varargin{4};
DomInfo      = varargin{5};
        
% [Fxi Feta] = force_oneform(r,c,FunctionType,Domain,DomInfo)
        
xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numRows+1)   ;   % This might be more advanced !!!!!

xiLR  = (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xi';
etaAB = (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*eta;

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

varargout{1} = Fxi;
varargout{2} = Feta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 2
        
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

F = zeros(N*N,1);
for i=1:N
    for j=1:N
        xib  = ( (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xibar )*ones(1,n);
        etab = ones(n,1)*( (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etabar' );
        
        [xx,yy,dxdxib,dxdetab,dydxib,dydetab] = coordinatemap(xib,etab,Domain,DomInfo);

        f = exact_solution(xx,yy,FunctionType,'force');

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
        F(ij,1) = sum(sum(GGw.*f.*J));
    end
end

varargout{1} = F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end