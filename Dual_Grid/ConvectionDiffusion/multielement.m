clear all
% close all
clc

global N m cc numRows numColumns

numRows    = 2;
numColumns = 2;
N = 2;
m = 1;

N2 = N*N;

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);
cells_in_element = N2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering unknowns                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localnr = reshape(1:N2,N,N);

globalnr = zeros(numColumns*N,numRows*N);
for j=1:numRows
    for i=1:numColumns
        globalnr((i-1)*N+(1:N),(j-1)*N+(1:N)) = ((i-1)+(j-1)*numColumns)*N2+localnr;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid option
cc = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL] = GLLnodes(N);    etaGLL = xiGLL;  % Gauss-Lobotto-Legendre
[xiG,wG]   = Gnodes(N);        etaG   = xiG;    % Gauss
xiEG     = [-1 xiG 1];         etaEG  = xiEG;   % Extended Gauss

XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';
XiGLLG  = xiGLL'*ones(1,N);    EtaGLLG   = ones(N+1,1)*etaG;
XiGGLL  = xiG'*ones(1,N+1);    EtaGGLL   = ones(N,1)*etaGLL;
XiEGEG  = xiEG'*ones(1,N+2);   EtaEGEG   = XiEGEG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numRows+1)   ;   % This might be more advanced !!!!!

% Global mesh
XGLLGLL = zeros(N*numColumns+1,N*numRows+1);
YGLLGLL = zeros(N*numColumns+1,N*numRows+1);
XGG   = zeros(N*numColumns  ,N*numRows  );
YGG   = zeros(N*numColumns  ,N*numRows  );
ind1 = 2*N*(N+1)*numColumns*numRows;
ind2 = 4*N^2*(N+1)^2*numColumns*numRows;
B = spalloc(ind1,ind1,ind2);
ind1 = (N2+3*N)*numColumns*numRows;
ind2 = (N2+3*N)^2*numColumns*numRows;
C = spalloc(ind1,ind1,ind2);
ind = numColumns*numRows*(N2+2*N)+(numRows+numColumns)*N;
f = zeros(ind,1);
for r=1:numRows
    for c=1:numColumns

        XibGLLGLL  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiGLL )'*ones(1,N+1);
        EtabGLLGLL = ones(N+1,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*etaGLL );

        XibGG  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiG )'*ones(1,N);
        EtabGG = ones(N,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*etaG );

        % gridtype = 'sinecurve' % This might be more advanced !!!!!
        XeGLLGLL = XibGLLGLL + cc*sin(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
        YeGLLGLL = EtabGLLGLL+ cc*sin(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);

        XeGG = XibGG + cc*sin(pi*XibGG).*sin(pi*EtabGG);
        YeGG = EtabGG+ cc*sin(pi*XibGG).*sin(pi*EtabGG);

        if c<numColumns && r<numRows
            XGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N)) = XeGLLGLL(1:N,1:N);
            YGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N)) = YeGLLGLL(1:N,1:N);
        elseif c<numColumns && r==numRows
            XGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = XeGLLGLL(1:N,1:N+1);
            YGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = YeGLLGLL(1:N,1:N+1);
        elseif c==numColumns && r<numRows
            XGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = XeGLLGLL(1:N+1,1:N);
            YGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = YeGLLGLL(1:N+1,1:N);
        elseif c==numColumns && r==numRows
            XGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = XeGLLGLL;
            YGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = YeGLLGLL;
        end
        XGG((c-1)*N+(1:N),(r-1)*N+(1:N)) = XeGG;
        YGG((c-1)*N+(1:N),(r-1)*N+(1:N)) = YeGG;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [Be,Ce] = elementmatrix_method2(N,xiGLL,xiEG,wG,wGLL,XibGLLGLL,EtabGLLGLL);
        ind = 2*N*(N+1)*((r-1)*numColumns+c-1)+(1:2*N*(N+1));
        B(ind,ind) = Be;

        ind = (c-1+(r-1)*numColumns)*(N2+4*N)+(1:N2+4*N);
        C(ind,ind) = Ce;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        F = force(XibGLLGLL(:,1),EtabGLLGLL(1,:));
        i = (c-1)*N+(1:N);
        j = (r-1)*N+(1:N);
        ind = globalnr(i,j)+...
              +(r-1)*(2*numColumns+1)*N+...
              +(r==1)*(c>1)*(c-1)*N+(r>1)*numColumns*N+...
              +(c-1)*2*N+(c>1)*N;
        f(ind) = F;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for r = numRows:-1:1
%     for c = numColumns:-1:1
%         max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
%         max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
%         max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;
% 
%         if r>1
%             for i=1:N
%                 ind1 = max_in_lower_element-2*(i-1);
%                 ind2 = max_in_element-(2*i-1);
% 
%             end
%         end
%         if c>1
%             for i=1:N
%                 ind1 = max_in_left_element-2*N-2*(i-1);
%                 ind2 = max_in_element-2*N-(2*i-1);
%                 
%                 
%             end
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dpe,Gde] = topology(N);

% Divergence operator

Dp = kron(speye(numRows*numColumns),Dpe);

for r = numRows:-1:1
    for c = numColumns:-1:1
        max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
        max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
        max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;

        if r>1
            for i=1:N
                ind1 = max_in_lower_element-2*(i-1);
                ind2 = max_in_element-(2*i-1);
                Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
                Dp(ind2,:) = [];
                C(:,ind2) = [];
                C(ind2,:) = [];
            end            
        end
        if c>1
            for i=1:N
                ind1 = max_in_left_element-2*N-2*(i-1);
                ind2 = max_in_element-2*N-(2*i-1);
                Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
                Dp(ind2,:) = [];
                C(:,ind2) = [];
                C(ind2,:) = [];
            end
        end
    end
end

% Gradient operator
Gd = -Dp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct matrix                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = -Dp*inv(B)*Dp'; %clear B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Hconv>1 %&& method == 1
% remove bc columns

max_in_element = (numRows*numColumns)*N2+(numRows*(numColumns-1)+numColumns*(numRows-1))*N+... % internal points
                 2*N*(numColumns+numRows); % boundary points
for r = numRows:-1:1
    for c = numColumns:-1:1
        if r == numRows
            for i=N:-1:1
                ind = max_in_element-(N-i)*(1+(r==1));
                A( : , ind ) = [];
                A( ind , : ) = [];
                f(ind) = [];
            end
        end
        if r == 1
            for i=N:-1:1
                ind = max_in_element-(1+2*(N-i))-(r==numRows)*(i-1);
                A( : , ind ) = [];
                A( ind , : ) = [];
                f(ind) = [];
            end
        end
        if c == numColumns
            for i=N:-1:1
                ind = max_in_element-N-(N-i)-N*(r==1)-(N-i)*(c==1);
                A( : , ind ) = [];
                A( ind , : ) = [];
                f(ind) = [];
            end
        end
        if c == 1
            for i=N:-1:1
                ind = max_in_element-N-1-2*(N-i)-N*(r==1)-(i-1)*(c==numColumns);
                A( : , ind ) = [];
                A( ind , : ) = [];
                f(ind) = [];
            end
        end
        max_in_element = max_in_element-(N2+2*N)-N*((c==1)+(r==1));
    end
end
end

%%%%%%%%%%%%%%%%%%%%%

% fbc = zeros((numRows*(numColumns-1)+numColumns*(numRows-1))*N,1);
% 
% f = [ f ; fbc];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = A\f;
% keyboard
if method == 2
    
    phi_in = phi;
    phi = zeros(numRows*numColumns*(N2+2*N)+(numRows+numColumns)*N,1);
    
    inner = [];
    for r=1:numRows
        for c=1:numColumns
            max_prev_element = ((r-1)*numColumns+c-1)*(N2+2*N)+((r==1)*(c-1)+(r>1)*numColumns)*N;
            inner = [inner max_prev_element+(r-1+(c>1))*N+(1:N2)];
            if c<numColumns
                inner = [inner max_prev_element+(r-1+(c>1))*N+N2+((1+(c==1)):(1+(c==1)):(1+(c==1))*N)];
            end
            if r<numRows
                inner = [inner max_prev_element+(r-1+(c>1))*N+N2+N+(c==1)*N+((1+(r==1)):(1+(r==1)):(1+(r==1))*N)];
            end
        end
    end
    
    phi(inner) = phi_in;
    
    phi = C\phi;
    
    phi = phi(inner);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen_me
