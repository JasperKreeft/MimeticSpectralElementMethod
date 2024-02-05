clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
if ispc
    
elseif isunix
run '/media/My Passport/MSEM/MSEM_codes/Single_Grid_V4/Convection/Skew/Library_Convection/GetLibrary.m'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns nn
global xi
global h e w
global nr_0 nr_1 nr_2
global globalnr_0 globalnr_1v globalnr_1h globalnr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'SinDeformGrid';
DomInfo = 0.0;

plot_fig  = 1;
avimatlab = 0;
Tecplot   = 0;

N = 5;
H = 5;
T = 4;
dt = 0.04;

numRows    = H;
numColumns = H;
numElements = numRows*numColumns;

filename = ['One_lin_ME_N' num2str(N) '_dt' num2str(dt)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize figure and plot Mesh

if ispc && (plot_fig || avimatlab)
    fig = figure('windowstyle','docked');
elseif  plot_fig || avimatlab
    fig = figure('Units','normalized','Position',[0, 0, 0.5, 1]);
end
if plot_fig
    set(gcf,'visible','on')
else
    set(gcf,'visible','off')
end

if plot_fig || avimatlab
    subplot(1,2,1)
    meshplot
    axis([min(min(Mesh.X)) max(max(Mesh.X)) min(min(Mesh.Y)) max(max(Mesh.Y))])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

U0x = zeros(N*(N+1),numElements);
U0y = zeros(N*(N+1),numElements);
for c=1:numColumns
    for r=1:numRows
        rc = c+(r-1)*numColumns;
        [U0x(:,rc),U0y(:,rc)] = initialforms(1,r,c,'gaussian',Domain,DomInfo);
    end
end

[UU0x,UU0y] = reconstruct(1,U0x,U0y,hp,ep,Meshp);

if plot_fig || avimatlab
    
    for i=1:numElements

        Xp = reshape(Meshp.X(:,i),nn,nn);
        Yp = reshape(Meshp.Y(:,i),nn,nn);

        subplot(1,2,1)
        pcolor(Xp,Yp,reshape(UU0x(:,i),nn,nn))
        hold on
        shading interp
        axis equal
        axis([-1 1 -1 1])
        set(gca,'clim',[-1 1])

        subplot(1,2,2)
        pcolor(Xp,Yp,reshape(UU0y(:,i),nn,nn))
        hold on
        shading interp
        axis equal
        axis([-1 1 -1 1])
        set(gca,'clim',[-1 1])

    end
    subplot(1,2,1); hold off
    subplot(1,2,2); hold off
    pause(0.15)
end

if Tecplot
    dataXY = [ Mesh.X Mesh.Y ];
    data = [ dataXY reshape(UU0x+UU0y,[],1) ];
    name = strcat(filename,'_00');
    MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization

NGe = normalgrad(N);
De  = div(N);

E01xe = zeros((N+1)^2,N*(N+1));
for i=1:N
    E01xe(:,(i-1)*(N+1)+(1:N+1)) = kron(eye(N+1),e(i,:)');
end
E01ye = kron(e',eye(N+1));



E12xe = kron(eye(N),e');
E12ye = zeros(N*(N+1),N*N);
for i=1:N
    E12ye(1:N*(N+1),(i-1)*N+(1:N)) = kron(eye(N),e(i,:)');
end

NG = zeros(nr_1,nr_0);
D  = zeros(nr_2,nr_1);
M1 = zeros(nr_1);
vE01x = zeros(nr_0,nr_1);
vE01y = zeros(nr_0,nr_1);
vE12x = zeros(nr_1,nr_2);
vE12y = zeros(nr_1,nr_2);

for i=1:numElements

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
ind2 = globalnr_2(:,i);

% one-forms
M1e = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
% M1(ind1,ind1) = M1(ind1,ind1) + M1e;
M1(ind1,ind1) = M1e;

% normal gradient operator
NG(ind1,ind0) = NGe;

% Divergence operator
D(ind2,ind1) = De;

Axe = 1*H;
Aye = 1*H;

ind1_v = globalnr_1v(:,i);
% vE01y(ind0,ind1_v) = vE01y(ind0,ind1_v) + Aye.*E01ye;
vE01y(ind0,ind1_v) = Aye.*E01ye;
% vE12x(ind1_v,ind2) = vE12x(ind1_v,ind2) + Axe.*E12xe;
vE12x(ind1_v,ind2) = Axe.*E12xe;

ind1_h = globalnr_1h(:,i);
% vE01x(ind0,ind1_h) = vE01x(ind0,ind1_h) + Axe.*E01xe;
vE01x(ind0,ind1_h) = Axe.*E01xe;
% vE12y(ind1_h,ind2) = vE12y(ind1_h,ind2) + Aye.*E12ye;
vE12y(ind1_h,ind2) = Aye.*E12ye;


end

vE01 = -vE01y + vE01x;
vE12 =  vE12x + vE12y;

%%

% Periodic BC
ind1L = unique(reshape(globalnr_1v(1:(N+1):N*(N+1),1:numColumns:numElements),[],1));
ind1R = unique(reshape(globalnr_1v(N+1:(N+1):N*(N+1),numColumns:numColumns:numElements),[],1));
ind1B = unique(reshape(globalnr_1h(1:N+1:N*(N+1),1:numColumns),[],1));
ind1T = unique(reshape(globalnr_1h(N+1:N+1:N*(N+1),(numRows-1)*numColumns+1:numElements),[],1));

NG([ind1L;ind1B],:) = NG([ind1L;ind1B],:)+NG([ind1R;ind1T],:);
NG([ind1R;ind1T],:) = [];

D(:,[ind1L;ind1B]) = D(:,[ind1L;ind1B])+D(:,[ind1R;ind1T]);
D(:,[ind1R;ind1T]) = [];

vE01(:,[ind1L;ind1B]) = (vE01(:,[ind1L;ind1B])+vE01(:,[ind1R;ind1T]))/2; % central discretization
vE01(:,[ind1R;ind1T]) = [];

vE12([ind1L;ind1B],:) = ( vE12([ind1L;ind1B],:)+vE12([ind1R;ind1T],:) )/2; % central discretization
vE12([ind1R;ind1T],:) = [];

M1([ind1L ind1B],:) = M1([ind1L ind1B],:)+M1([ind1R;ind1T],:);
M1(:,[ind1L ind1B]) = M1(:,[ind1L ind1B])+M1(:,[ind1R;ind1T]);
M1([ind1R;ind1T],:) = [];
M1(:,[ind1R;ind1T]) = [];

ind0B = unique(reshape(globalnr_0(1:N+1,1:numColumns),[],1));
ind0T = unique(reshape(globalnr_0(N*(N+1)+(1:N+1),(numRows-1)*numColumns+1:numElements),[],1));
ind0L = unique(reshape(globalnr_0(1:N+1:(N+1)^2,1:numColumns:numElements),[],1)); ind0L(end) = [];
ind0R = unique(reshape(globalnr_0(N+1:N+1:(N+1)^2,numColumns:numColumns:numElements),[],1)); ind0R(end) = [];

NG(:,[ind0L;ind0B]) = NG(:,[ind0L;ind0B])+NG(:,[ind0R;ind0T]);
NG(:,[ind0R;ind0T]) = [];

ind_element_boundary = [1:N+1 N+2:N+1:N*(N+1) 2*(N+1):N+1:N*(N+1) N*(N+1)+1:(N+1)^2];
ind_all_el_bcs = unique(reshape(globalnr_0(ind_element_boundary,:),[],1));
vE01(ind_all_el_bcs,:) = vE01(ind_all_el_bcs,:)/2; % central discretization

ind = unique(globalnr_0([1 N+1 N*(N+1)+1 (N+1)^2],:));
vE01(ind,:) = vE01(ind,:)/2; % central discretization

vE01(ind0B,:) = (vE01(ind0B,:)+vE01(ind0T,:));
vE01(ind0L,:) = (vE01(ind0L,:)+vE01(ind0R,:));
vE01([ind0R;ind0T],:) = [];



% Periodic BC

% Vertical cell-interfaces
ind = 0;
for r=1:numRows
    for c=1:numColumns-1
        rc = c+(r-1)*numColumns;

        ind = ind + (c==1)*( (r>1)*(N*(N-1)+N*N) + (r==2)*N );
        
        ind1 = ind + ((N+(c==1)):(N+(c==1)):N*(N+(c==1)));

        vE12(ind1,:) = vE12(ind1,:)/2; % central discretization

        ind = ind+(N*N+N*(N-1+(r<numRows)))+(r==1)*N+(c==1)*N;

    end
end


% Horizontal cell-interfaces
ind = 0;
for r=1:numRows-1
    for c=1:numColumns

        ind = ind + N*(N-1+(c<numColumns))+(c==1)*N;
        
        ind1 = ind + ((N+(r==1)):(N+(r==1)):N*(N+(r==1)));

        vE12(ind1,:) = vE12(ind1,:)/2; % central discretization

        ind = ind+(N*(N-1+(r<numRows)))+(r==1)*N;

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NG = (NG>0)-(NG<0);
% ans = full(NG*vE01);
% break
A = M1*(vE12*D+NG*vE01);
A = -A';

% K1 = (A-A')/2; % Skew-symmetric
K1 = A;      % Standard
K2 = M1;

K1 = sparse(K1);
K2 = sparse(K2);

U = zeros(nr_1,1);
U(globalnr_1v) = U0x;
U(globalnr_1h) = U0y;
U([ind1R;ind1T]) = [];

t = 0; j=0;
while t<T
    j = j+1
    t = t+dt;
    
    U = timemarching(K2,K1,dt,U,'ESDIRK3');
    
ind = 1:nr_1; ind([ind1R;ind1T]) = [];

    UU = zeros(nr_1,1);
    UU(ind) = U;
    UU(ind1R) = UU(ind1L);
    UU(ind1T) = UU(ind1B);
    Ux = UU(globalnr_1v);
    Uy = UU(globalnr_1h);
    
    [UUx,UUy] = reconstruct(1,Ux,Uy,hp,ep,Meshp);
    
    if plot_fig || avimatlab
        for i=1:numElements

            Xp = reshape(Meshp.X(:,i),nn,nn);
            Yp = reshape(Meshp.Y(:,i),nn,nn);

            subplot(1,2,1)
            pcolor(Xp,Yp,reshape(UUx(:,i),nn,nn))
            hold on
            shading interp
            axis equal
            axis([-1 1 -1 1])
            set(gca,'clim',[-1 1])

            subplot(1,2,2)
            pcolor(Xp,Yp,reshape(UUy(:,i),nn,nn))
            hold on
            shading interp
            axis equal
            axis([-1 1 -1 1])
            set(gca,'clim',[-1 1])

        end
        subplot(1,2,1); hold off
        subplot(1,2,2); hold off
        pause(0.05)
    end
    
    if avimatlab
        avimovie(filename,fig,t==dt,t>=T);
    end

    if Tecplot
        data = [ dataXY reshape(UUx+UUy,[],1) ];
        if j<10
            name = strcat(filename,'_0',num2str(j));
        else
            name = strcat(filename,'_',num2str(j));
        end
        MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
if ispc
    
elseif isunix
run '/media/My Passport/MSEM/MSEM_codes/Single_Grid_V4/Convection/Skew/Library_Convection/GetLibrary.m'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%