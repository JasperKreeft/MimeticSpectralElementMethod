clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test 1: cosine hill
% test 2: full cosine hill
% test 3: step function
% test 4: cosine function

test = 1;

N = 4;
N2 = N*N;

numColumns = 16;
numSlabs = numColumns;

a = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering unknowns                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localnr = reshape(1:N2,N,N);

globalnr = zeros(numColumns*(N+1)+1,N+2);
for c=1:numColumns
    globalnr((c-1)*(N+1)+1+(1:N),2:N+1) = (c-1)*(N2+3*N)+N*(c>1)+localnr;
    if c==1
        globalnr(1,2:N+1)   = N2+(1:2:2*N);
        globalnr(N+2,2:N+1) = N2+(2:2:2*N);
    else
        globalnr(c*(N+1)+1,2:N+1) = (c-1)*(N2+2*N)+c*N+N2+(1:N);
    end
    globalnr((c-1)*(N+1)+1+(1:N),1)   = (c-1)*(N2+3*N)+(N2+2*N)+(1:2:2*N);
    globalnr((c-1)*(N+1)+1+(1:N),N+2) = (c-1)*(N2+3*N)+(N2+2*N)+(2:2:2*N);
end

% disp('  '); disp('globalnr = '); disp('  ');
% disp(num2str(flipud(globalnr')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiG,wG]   = Gnodes(N);        etaG   = xiG;    % Gauss
xiEG     = [-1 xiG 1];         etaEG  = xiEG;   % Extended Gauss

XiEGEG  = xiEG'*ones(1,N+2);   EtaEGEG   = XiEGEG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numSlabs+1);     % This might be more advanced !!!!!

% Global mesh
XGG   = zeros(N*numColumns  ,N );
YGG   = zeros(N*numSlabs  ,N );
XEGEG = zeros((N+1)*numColumns+1,N+2);
YEGEG = zeros((N+1)*numSlabs+1,N+2);

for t=1:numSlabs
    for c=1:numColumns

        XibGG  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiG )'*ones(1,N);
        EtabGG = ones(N,1)*( (etabAB(t+1)+etabAB(t))/2+(etabAB(t+1)-etabAB(t))/2*etaG );

        XeGG = XibGG;
        YeGG = EtabGG;

        XGG((c-1)*N+(1:N),(t-1)*N+(1:N)) = XeGG;
        YGG((c-1)*N+(1:N),(t-1)*N+(1:N)) = YeGG;

        XibEGEG  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiEG )'*ones(1,N+2);
        EtabEGEG = ones(N+2,1)*( (etabAB(t+1)+etabAB(t))/2+(etabAB(t+1)-etabAB(t))/2*etaEG );

        XeEGEG = XibEGEG;
        YeEGEG = EtabEGEG;

        XEGEG((c-1)*(N+1)+(1:N+2),(t-1)*(N+1)+(1:N+2)) = XeEGEG;
        YEGEG((c-1)*(N+1)+(1:N+2),(t-1)*(N+1)+(1:N+2)) = YeEGEG;

    end
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C1 = spalloc(nr_of_nodes,2*(N2+2*N),2*(N2+2*N));

for j=1:N
ind1 = (j-1)*N+(1:N);
ind2 = (j-1)*(N+2)+1+(1:N);
C1(ind1,ind2) = speye(N);
end

ind1 = N2+(1:2:2*N);
ind2 = 1:N+2:N*(N+2);
C1(ind1,ind2) = speye(N);

ind1 = N2+(2:2:2*N);
ind2 = (1:N)*(N+2);
C1(ind1,ind2) = speye(N);

for i=1:N
ind1 = (1:N:N2)+(i-1);
ind2 = N*(N+2)+1+(1:N)+(i-1)*(N+2);
C1(ind1,ind2) = speye(N);
end

ind1 = N2+2*N+(1:2:2*N);
ind2 = N*(N+2)+(1:N+2:N*(N+2));
C1(ind1,ind2) = speye(N);

ind1 = N2+2*N+1+(1:2:2*N);
ind2 = N*(N+2)+(1:N)*(N+2);
C1(ind1,ind2) = speye(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IT = C1*IT;

IT = kron(speye(numColumns),IT);

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);
cells_in_element = N2;

for c=numColumns:-1:2

    max_in_left_element  = (c-1)*nodes_in_element;

    ind1 = max_in_left_element+N2+(1:2:2*N);
    ind2 = max_in_left_element-4*N+(2:2:2*N);
    IT(ind2,:) = IT(ind2,:)+IT(ind1,:);
    IT(ind1,:) = [];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dpe,Gde] = topology(N);

% Divergence operator

Dp = kron(speye(numColumns),Dpe);

for c = numColumns:-1:1
    max_in_element       = c*nodes_in_element;
    max_in_left_element  = (c-1)*nodes_in_element;

    if c>1
        for i=1:N
            ind1 = max_in_left_element-2*N-2*(i-1);
            ind2 = max_in_element-2*N-(2*i-1);
            Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
            Dp(ind2,:) = [];
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
end

B = A;
ind_leftbc  = globalnr(1,:) ; ind_leftbc(ind_leftbc==0)   = [];
ind_lowerbc = globalnr(:,1)'; ind_lowerbc(ind_lowerbc==0) = [];
B([ind_leftbc ind_lowerbc],:) = [];
B_bcic = B(:,[ind_leftbc ind_lowerbc]);
B(:,[ind_leftbc ind_lowerbc]) = [];


indG = [];
for c=1:numColumns
    indG = [indG (c-1)*(N+1)+1+(1:N)];
end

for t=1:numSlabs

phi_bc = zeros(N,1);
if test==4
    % cosine function
    phi_bc = 1/4+1/4*cos(-2*pi*a*(YGG(1,(t-1)*N+(1:N))'+1));
end

f=-B_bcic*[phi_bc ; phi_ic];

phi = B\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PHI = zeros(numColumns*(N+1)+1,N+2);
for c=1:numColumns
    ind = (c-1)*(N2+2*N)+(1:N2);
    PHI((c-1)*(N+1)+1+(1:N),2:N+1) = reshape(phi(ind),N,N);

    ind = (c-1)*(N2+2*N)+N2+(1:N);
    PHI(c*(N+1)+1,2:N+1) = reshape(phi(ind),1,N);

    ind = (c-1)*(N2+2*N)+(N2+N)+(1:N);
    PHI(1+(c-1)*(N+1)+(1:N),N+2) = reshape(phi(ind),N,1);
end

PHI(indG,1) = phi_ic;     % Lower bc

PHI(1,2:N+1) = phi_bc';  % Left bc

phi_ic = PHI(indG,N+2);
% subplot(2,1,1); plot(XGG(:,1),phi_ic,'-xr'); ylim([-0.2 0.7])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = 50;
% [xx,wg] = Gnodes(nn); yy=xx;
xx = linspace(-1,1,nn); yy=xx;
h_p = LagrangeVal(xx,nn,2);
X = xx'*ones(1,nn);
Y = ones(nn,1)*yy;

hG   = LagrangeVal(xx,N,2);
hEG  = LagrangeVal(xx,N,3);
hGEG = (hG+hEG(2:N+1,:))/2;

pphi = zeros(nn,nn,numColumns);
pphiG = zeros(nn,nn,numColumns);

% subplot(2,1,2)
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 0.6])
% set(gca,'clim',[-.1 1.1])
axis equal
axis([-1 1 -1 1])
hold on

Yrc = ones(nn,1)*((etabAB(t+1)+etabAB(t))/2+(etabAB(t+1)-etabAB(t))/2*yy);
yrc = (etabAB(t+1)+etabAB(t))/2+(etabAB(t+1)-etabAB(t))/2*etaEG;

for c=1:numColumns

    rc = c;

    Xrc = ((xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xx)'*ones(1,nn);
    xrc = (xibLR(c+1)+xibLR(c))/2+(xibLR(c+1)-xibLR(c))/2*xiEG;

    XXrc = Xrc;
    YYrc = Yrc;

    xxrc = xrc'*ones(1,N+2);
    yyrc = ones(N+2,1)*yrc ;

    pphi(:,:,rc) = hEG(2:N+1,:)'*PHI((c-1)*(N+1)+1+(1:N),2:N+1)*hEG(2:N+1,:)+... % inner part
                   hEG(1,:)'*PHI((c-1)*(N+1)+1,2:N+1)*hGEG+...      % left side
                   hEG(N+2,:)'*PHI(c*(N+1)+1,2:N+1)*hGEG+...        % right side
                   hGEG'*PHI((c-1)*(N+1)+1+(1:N),1)*hEG(1,:)+...    % lower side
                   hGEG'*PHI((c-1)*(N+1)+1+(1:N),N+2)*hEG(N+2,:);   % upper side
               
    pphiG(:,:,rc) = hG'*PHI((c-1)*(N+1)+1+(1:N),2:N+1)*hG;

    surf(XXrc,YYrc,pphi(:,:,rc))%,'FaceColor','interp','EdgeColor','none','FaceLighting','phong')
%     surf(XXrc,YYrc,pphiG(:,:,rc))

end
shading interp
% view(-22,50)
% camlight left
% view([0 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pause(0.5)
end