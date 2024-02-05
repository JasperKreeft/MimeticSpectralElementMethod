clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
if ispc
% run 'O:\MSEM\MSEM_codes\Single_Grid_V5\Convection\Skew\Library_Convection\GetLibrary.m'
elseif isunix
run '/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Convection/Skew/Library_Convection/GetLibrary.m'
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Topology
% 
% numbering('square');
% 
% NGe = normalgrad(N);
% De  = div(N);
% 
% NG = zeros(nr_1,nr_0);
% D  = zeros(nr_2,nr_1);
% for i=1:numElements
% 
% ind0 = globalnr_0(:,i);
% ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
% ind2 = globalnr_2(:,i);
% 
% % normal gradient operator
% NG(ind1,ind0) = NGe;
% 
% % Divergence operator
% D(ind2,ind1) = De;
% 
% end
% 
% % Periodic BC
% ind1L = unique(reshape(globalnr_1v(1:(N+1):N*(N+1),1:numColumns:numElements),[],1));
% ind1R = unique(reshape(globalnr_1v(N+1:(N+1):N*(N+1),numColumns:numColumns:numElements),[],1));
% ind1B = unique(reshape(globalnr_1h(1:N+1:N*(N+1),1:numColumns),[],1));
% ind1T = unique(reshape(globalnr_1h(N+1:N+1:N*(N+1),(numRows-1)*numColumns+1:numElements),[],1));
% 
% NG([ind1L;ind1B],:) = NG([ind1L;ind1B],:)+NG([ind1R;ind1T],:);
% NG([ind1R;ind1T],:) = [];
% 
% D(:,[ind1L;ind1B]) = D(:,[ind1L;ind1B])+D(:,[ind1R;ind1T]);
% D(:,[ind1R;ind1T]) = [];
% 
% ind0B = unique(reshape(globalnr_0(1:N+1,1:numColumns),[],1));
% ind0T = unique(reshape(globalnr_0(N*(N+1)+(1:N+1),(numRows-1)*numColumns+1:numElements),[],1));
% ind0L = unique(reshape(globalnr_0(1:N+1:(N+1)^2,1:numColumns:numElements),[],1)); ind0L(end) = [];
% ind0R = unique(reshape(globalnr_0(N+1:N+1:(N+1)^2,numColumns:numColumns:numElements),[],1)); ind0R(end) = [];
% 
% NG(:,[ind0L;ind0B]) = NG(:,[ind0L;ind0B])+NG(:,[ind0R;ind0T]);
% NG(:,[ind0R;ind0T]) = [];
% 
% NG = (NG>0)-(NG<0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square_periodic');

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
% l=0; m=0;
% for j=ind1'
%     l=l+1;
%     for k=ind0'
%         m=m+1;
%         NG(j,k) = NG(j,k) + NGe(l,m);
%     end
%     m=0;
% end

% for l=1:length(unique(ind1))
%     j=ind1(l);
%     for k=ind0'
%         m=m+1;
%         NG(j,k) = NG(j,k) + NGe(l,m);
%     end
%     m=0;
% end






% Divergence operator
% D(ind2,ind1) = De;
k = 0;
for j=ind1'
    k=k+1;
D(ind2,j) = D(ind2,j) + De(:,k);
end

Axe = 1*H;
Aye = 1*H;

ind1_v = globalnr_1v(:,i);
vE01y(ind0,ind1_v) = Aye.*E01ye;
vE12x(ind1_v,ind2) = Axe.*E12xe;

ind1_h = globalnr_1h(:,i);
vE01x(ind0,ind1_h) = Axe.*E01xe;
vE12y(ind1_h,ind2) = Aye.*E12ye;

end

vE01 = -vE01y + vE01x;
vE12 =  vE12x + vE12y;



% NG = [  1  0 -1  0
%         0  1  0 -1
%        -1  0  1  0
%         0 -1  0  1
%        -1  1  0  0
%         0  0 -1  1
%         1 -1  0  0
%         0  0  1 -1  ];
% 
% 
% M1 = [  0.777777777777778	0	-0.111111111111111	0	0	0	0	0
%         0	1.55555555555556	0	-0.222222222222222	0	0	0	0
%        -0.111111111111111	0	0.777777777777778	0	0	0	0	0
%         0	-0.222222222222222	0	1.55555555555556	0	0	0	0
%         0	0	0	0	0.777777777777778	0	-0.111111111111111	0
%         0	0	0	0	0	1.55555555555556	0	-0.222222222222222
%         0	0	0	0	-0.111111111111111	0	0.777777777777778	0
%         0	0	0	0	0	-0.222222222222222	0	1.55555555555556   ];
% 
% 
% 
% vE12 = [    0.5000    0.5000         0         0
%             0.5000    0.5000         0         0
%                  0         0    0.5000    0.5000
%                  0         0    0.5000    0.5000
%             0.5000         0    0.5000         0
%             0.5000         0    0.5000         0
%                  0    0.5000         0    0.5000
%                  0    0.5000         0    0.5000 ];
%     
% vE01 = [   -0.2500         0   -0.2500         0    0.2500         0    0.2500         0
%                  0   -0.5000         0   -0.5000    0.2500         0    0.2500         0
%            -0.2500         0   -0.2500         0         0    0.5000         0    0.5000
%                  0   -0.5000         0   -0.5000         0    0.5000         0    0.5000 ];

% close; clc
% full(NG)
% keyboard



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = M1*(vE12*D+NG*vE01);
A = -A';

K1 = (A-A')/2; % Skew-symmetric   % Werkt niet!!!
% K1 = A;      % Standard
% clear M1
% load M1_N23H1.mat
K2 = M1;

K1 = sparse(K1);
K2 = sparse(K2);

U = zeros(nr_1,1);
U(globalnr_1v) = U0x;
U(globalnr_1h) = U0y;

t = 0; j=0;
while t<T
    j = j+1
    t = t+dt;
    
    U = timemarching(K2,K1,dt,U,'ESDIRK3');
    
    Ux = U(globalnr_1v);
    Uy = U(globalnr_1h);
    
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
run '/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Convection/Skew/Library_Convection/GetLibrary.m'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%