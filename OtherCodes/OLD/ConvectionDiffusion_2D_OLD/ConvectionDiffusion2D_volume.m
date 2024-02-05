clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N N2 numRows numColumns RC RCe
global nodes_in_element edges_in_element cells_in_element

numColumns = 5;
numRows    = 5;

N = 5;
N2 = N*N;

global m
m=1;

Re = 1;
velocity = [ 0 ; 0 ];

global BC

%      le ri lo up
BC = [ 0  0  0  1      % q
       1  1  1  0 ];   % phi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);
cells_in_element = N2;

RC  = numColumns*numRows;
RCe = RC*edges_in_element;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering of 0-cells                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XGLLGLL,YGLLGLL,XGG,YGG,QGLLGLL,JGLLGLL,JGG] = gridgenerator();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global xiGLL xiEG
global h h_w hw_w e e_w ew_w

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);
[hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);
e    = EdgeVal(dhdxi );
e_w  = EdgeVal(dhdxiw);
ew_w = EdgeVal(dhwdxiw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T1,T2] = reorder();
    
    ind1 = nodes_in_element*RC;
    ind2 = edges_in_element*RC;
    
    ind = edges_in_element^2*RC;
S1 = spalloc(ind2,ind2,ind);

    ind = (nodes_in_element)^2*RC; % might be more accurate
S0H02 = spalloc(ind1,ind1,ind);

    ind = RC*2*N2*(N+1);
X12  = spalloc(ind2,ind1,ind);

for r=1:numRows
    for c=1:numColumns
        rc = (c+(r-1)*numColumns);

        % Diffusive part
        [S1e,S0e,W0] = ElementSupportOperatorMethod(QGLLGLL(:,3*(rc-1)+(1:3)),JGLLGLL(:,rc));
        [H02e]       = ElementHodgeMatrix(JGG(1:N2,rc));

        S0e(1:N2,1:N2) = S0e(1:N2,1:N2)*H02e;

        % Convective part
        X12e = ElementConvectionVolumeMatrix(velocity);

        % Assembly
            ind1 = (rc-1)*nodes_in_element+(1:nodes_in_element);
            ind2 = (rc-1)*edges_in_element+(1:edges_in_element);
            ind3 = (rc-1)*nodes_in_element+(1:cells_in_element);
        S1(ind2,ind2)    = S1e;
        S0H02(ind1,ind1) = S0e;
        X12(ind2,ind3)   = X12e;
    end
end

W0 = kron(speye(RC),W0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dpe = topology(N);

% Divergence operator

Dp = kron(speye(RC),Dpe);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = numRows:-1:1
    for c = numColumns:-1:1
        max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
        max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
        max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;

        if r>1
            ind1 = max_in_lower_element-2*((1:N)-1);
            ind2 = max_in_element-(2*(1:N)-1);
            S0H02(:,ind2) = [];
            S0H02(ind2,:) = [];
            W0(:,ind2) = [];
            W0(ind2,:) = [];
            Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
            Dp(ind2,:) = [];
            X12(:,ind1) = X12(:,ind1) + X12(:,ind2);
            X12(:,ind2) = [];
        end
        if c>1
            ind1 = max_in_left_element-2*N-2*((1:N)-1);
            ind2 = max_in_element-2*N-(2*(1:N)-1);
            S0H02(:,ind2) = [];
            S0H02(ind2,:) = [];
            W0(:,ind2) = [];
            W0(ind2,:) = [];
            Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
            Dp(ind2,:) = [];
            X12(:,ind1) = X12(:,ind1) + X12(:,ind2);
            X12(:,ind2) = [];
        end
    end
end

% Gradient operator
Gd = -Dp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct matrix                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A = [ S1/Re    Dp'*S0H02
%       Dp zeros(RC*(N2+2*N)+(numRows+numColumns)*N) ];

% A = [ -eye(size(Dp,2)) X12
%         Dp     zeros(RC*(N2+2*N)+(numRows+numColumns)*N) ];

A = [ S1/Re       Dp'*S0H02+X12
       Dp   zeros(RC*(N2+2*N)+(numRows+numColumns)*N) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind1 = RC*edges_in_element; % number of edges
ind2 = RC*(N2+2*N)+(numRows+numColumns)*N; % number of points
f1 = zeros(ind1,1);
f2 = zeros(ind2,1);
for r=1:numRows
    for c=1:numColumns

        F = force(r,c);
        i = 1+(c-1)*(N+1)+(1:N);
        j = 1+(r-1)*(N+1)+(1:N);
        f2(globalnr_0(i,j)) = F;

    end
end

f = [ f1
      f2 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[qphi_bc,q_bc,phi_bc,ind_bc,sindq]=boundaryconditions(globalnr_0,globalnr_1h,globalnr_1v,XGLLGLL,YGLLGLL,XGG,YGG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_bc = A(:,ind_bc);
A(ind_bc,:) = [];
A(:,ind_bc) = [];

f_bc = -A_bc*qphi_bc;
f = f+f_bc;
f(ind_bc) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qphi = A\f;

ind = RC*edges_in_element-sindq;
q_in = qphi(1:ind);

p = qphi(ind+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% multi_element_volume_postprocessen

etaGLL = xiGLL;
dxideta = diff(xiGLL)'*diff(etaGLL);

ind = 0;
for r=1:numRows
for c=1:numColumns
    rc = c+(r-1)*numColumns;
    
    P = reshape(p(ind+(1:N2)),N,N);

    ind = ind+N2+N*((c<numColumns)+(r<numRows));

    subplot(1,2,1)
    for i=1:N
        for j=1:N
            ii = (c-1)*N+i;
            jj = (r-1)*N+j;
    surf([XGLLGLL(ii:ii+1,1)' ; XGLLGLL(ii:ii+1,1)'],[YGLLGLL(1,jj) YGLLGLL(1,jj) ; YGLLGLL(1,jj+1) YGLLGLL(1,jj+1)],P(i,j)/dxideta(i,j)*ones(2))
    hold on
    view([0 0 1])
        end
    end

nn = 50;
xx = linspace(-1,1,nn);
XX = linspace(-1+2/numColumns*(c-1),-1+2/numColumns*c,nn)'*ones(1,nn);
yy = linspace(-1,1,nn);
YY = ones(nn,1)*linspace(-1+2/numRows*(r-1),-1+2/numRows*r,nn);

XiGLL  = xiGLL'*ones(1,N+1);
EtaGLL = ones(N+1,1)*etaGLL;

[hh,dhhdxx] = LagrangeVal(xx,N,1);
ee = EdgeVal(dhhdxx);
pp = ee'*P*ee;

subplot(1,2,2)
surf(XX,YY,pp)
hold on
view([0 0 1])
end
end

subplot(1,2,1)
xlabel('x')
ylabel('y')
view([0 0 1])
axis square

subplot(1,2,2)
shading interp
xlabel('x')
ylabel('t')
% set(gca,'clim',[-1.1 1.1])
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen_q
break
phi_in = p;
postprocessen_me

errorL2

errorL2_interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%