clear all
close all
clc

global N N2 numRows numColumns
global nodes_in_element edges_in_element cells_in_element

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test 1: cosine hill
% test 2: full cosine hill
% test 3: step function
% test 4: cosine function

test = 4;

N = 8;
N2 = N*N;

numColumns = 5;
numRows = numColumns;

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);
cells_in_element = N2;

a = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering unknowns                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr = numbering();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiG,wG] = Gnodes(N);        etaG  = xiG;    % Gauss
xiEG     = [-1 xiG 1];       etaEG = xiEG;   % Extended Gauss

XiEGEG  = xiEG'*ones(1,N+2);   EtaEGEG   = XiEGEG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numRows+1)   ;   % This might be more advanced !!!!!

% Global mesh
XGG   = zeros(N*numColumns  ,N*numRows  );
YGG   = zeros(N*numColumns  ,N*numRows  );
XEGEG = zeros((N+1)*numColumns+1,(N+1)*numRows+1);
YEGEG = zeros((N+1)*numColumns+1,(N+1)*numRows+1);

for r=1:numRows
    for c=1:numColumns

        XibGG  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiG )'*ones(1,N);
        EtabGG = ones(N,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*etaG );

        XeGG = XibGG;
        YeGG = EtabGG;

        XGG((c-1)*N+(1:N),(r-1)*N+(1:N)) = XeGG;
        YGG((c-1)*N+(1:N),(r-1)*N+(1:N)) = YeGG;

        XibEGEG  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiEG )'*ones(1,N+2);
        EtabEGEG = ones(N+2,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*etaEG );

        XeEGEG = XibEGEG;
        YeEGEG = EtabEGEG;

        XEGEG((c-1)*(N+1)+(1:N+2),(r-1)*(N+1)+(1:N+2)) = XeEGEG;
        YEGEG((c-1)*(N+1)+(1:N+2),(r-1)*(N+1)+(1:N+2)) = YeEGEG;
        
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % forcing function is empty. Boundary part is added later

    end
end
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        
        [hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);
        ew_w           = EdgeVal(dhwdxiw);

        IT = kron(speye(N),ew_w');

        nr_of_nodes = N2+4*N;
        nr_of_lines = 2*N*(N+1);

        IT = [ a*IT zeros(N2+2*N,N*(N+1))
               zeros(N2+2*N,N*(N+1)) IT ];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [C1,C2] = reorder();

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        IT = C1*IT;
        
        IT = kron(speye(numColumns*numRows),IT);

        IT = C2*IT; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dpe = topology(N);

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
            end            
        end
        if c>1
            for i=1:N
                ind1 = max_in_left_element-2*N-2*(i-1);
                ind2 = max_in_element-2*N-(2*i-1);
                Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
                Dp(ind2,:) = [];
            end
        end
    end
end

% Gradient operator
Gd = -Dp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct matrix                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = IT*Gd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary / initial conditions                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_ic = zeros(numColumns*N,1);
phi_bc = zeros(numRows*N,1);
if test==1
    % cosine hill
    for i=1:numColumns*N
        if XGG(i,1)<=0
            phi_ic(i,1) = 1/4-1/4*cos(2*pi*XGG(i,1));
        else
            phi_ic(i,1) = 0;
        end
    end
elseif test==2
    % full cosine hill
    phi_ic = 1/4+1/4*cos(pi*XGG(:,1));
elseif test==3
    % step-function
    xstep = -0.5;
    phi_ic = (XGG(:,1)>xstep);
elseif test==4
    % cosine function
    phi_ic = 1/4+1/4*cos(2*pi*XGG(:,1));
    phi_bc = 1/4+1/4*cos(-2*pi*a*(YGG(1,:)'+1));
end


plot(XGG(:,1),phi_ic,'-x')
hold on
plot(YGG(1,:),phi_bc,'-or')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = A;
ind_leftbc  = globalnr(1,:) ; ind_leftbc(ind_leftbc==0)   = [];
ind_lowerbc = globalnr(:,1)'; ind_lowerbc(ind_lowerbc==0) = [];
B([ind_leftbc ind_lowerbc],:) = [];
B_bcic = B(:,[ind_leftbc ind_lowerbc]);
B(:,[ind_leftbc ind_lowerbc]) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve system                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=-B_bcic*[phi_bc ; phi_ic];

phi = B\f;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PHI = zeros(numColumns*(N+1)+1,numRows*(N+1)+1);
for r=1:numRows
    for c=1:numColumns
        ind = (r-1)*(numColumns*(N2+2*N))+(c-1)*(N2+2*N)+(1:N2);
        PHI((c-1)*(N+1)+1+(1:N),(r-1)*(N+1)+1+(1:N)) = reshape(phi(ind),N,N);

            ind = (r-1)*(numColumns*(N2+2*N))+(c-1)*(N2+2*N)+N2+(1:N);
            PHI(c*(N+1)+1,(1:N)+1+(r-1)*(N+1)) = reshape(phi(ind),1,N);

            ind = (r-1)*(numColumns*(N2+2*N))+(c-1)*(N2+2*N)+(N2+N)+(1:N);
            PHI(1+(c-1)*(N+1)+(1:N),1+r*(N+1) ) = reshape(phi(ind),N,1);
    end
end

ind = [];
for c=1:numColumns
    ind = [ind (c-1)*(N+1)+1+(1:N)];
end
PHI(ind,1) = phi_ic;              % Lower bc

ind = [];
for r=1:numRows
    ind = [ind (r-1)*(N+1)+1+(1:N)];
end
PHI(1,ind) = phi_bc';%zeros(1,numRows*N);  % Left bc

% flipud(PHI')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = 5;
xx = linspace(-1,1,nn); yy=xx;
h_p = LagrangeVal(xx,nn,2);
X = xx'*ones(1,nn);
Y = ones(nn,1)*yy;

hG   = LagrangeVal(xx,N,2);
hEG  = LagrangeVal(xx,N,3);
hGEG = (hG+hEG(2:N+1,:))/2;

pphi = zeros(nn,nn,numColumns*numRows);
pphiG = zeros(nn,nn,numColumns*numRows);

figure

for r=1:numRows
    
    Yrc = ones(nn,1)*((etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*yy);
    yrc = (etabAB(r+1)+etabAB(r))/2+(etabAB(r+1)-etabAB(r))/2*etaEG;
    
    for c=1:numColumns
        
        rc = c+(r-1)*numColumns;
        
        Xrc = ((xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xx)'*ones(1,nn);
        xrc = (xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xiEG;

        XXrc = Xrc;
        YYrc = Yrc;
        
        xxrc = xrc'*ones(1,N+2);
        yyrc = ones(N+2,1)*yrc ;
        
        pphi(:,:,rc) = hEG(2:N+1,:)'*PHI((c-1)*(N+1)+1+(1:N),(r-1)*(N+1)+1+(1:N))*hEG(2:N+1,:)+... % inner part
                       hEG(1,:)'*PHI((c-1)*(N+1)+1,(r-1)*(N+1)+1+(1:N))*hGEG+...    % left side
                       hEG(N+2,:)'*PHI(c*(N+1)+1,(r-1)*(N+1)+1+(1:N))*hGEG+...      % right side
                       hGEG'*PHI((c-1)*(N+1)+1+(1:N),(r-1)*(N+1)+1)*hEG(1,:)+...    % lower side
                       hGEG'*PHI((c-1)*(N+1)+1+(1:N),r*(N+1)+1)*hEG(N+2,:);         % upper side
                   
        pphiG(:,:,rc) = hG'*PHI((c-1)*(N+1)+1+(1:N),(r-1)*(N+1)+1+(1:N))*hG;

        surf(XXrc,YYrc,pphi(:,:,rc))%,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
%         surf(XXrc,YYrc,pphiG(:,:,rc))%,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
        hold on
        
    end
end

shading interp
xlabel('x')
ylabel('y')
colorbar
set(gca,'clim',[-.1 0.6])
% set(gca,'clim',[-.1 1.1])
axis equal
axis([-1 1 -1 1])
% view(-22,50)
% camlight left

view([0 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
