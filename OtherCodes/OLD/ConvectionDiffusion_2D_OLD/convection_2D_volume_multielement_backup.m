clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%

global N N2 numRows numColumns RC
global nodes_in_element edges_in_element cells_in_element

numColumns = 2;
numRows    = 1;
RC = numColumns*numRows;

N = 2;

N2 = N*N;

a = 2;
b = 1;
velocity = [ a ; b ];

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);
cells_in_element = N2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering of 0-cells                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dpe = topology(N);

% Divergence operator

Dp = kron(speye(RC),Dpe);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XGLLGLL,YGLLGLL,XGG,YGG,QGLLGLL,JGLLGLL,JGG] = gridgenerator();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global xiGLL
global h e

[h,dhdxi] = LagrangeVal(xiGLL,N,1);
e         = EdgeVal(dhdxi );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X12  = spalloc(RC*edges_in_element,RC*cells_in_element,RC*2*N2*(N+1));

for r=1:numRows
    for c=1:numColumns
        rc = (c+(r-1)*numColumns);

        % Convective part
        X12e = ElementConvectionVolumeMatrix(velocity);

        % Assembly
            ind1 = (rc-1)*edges_in_element+(1:edges_in_element);
            ind2 = (rc-1)*cells_in_element+(1:cells_in_element);
        X12(ind1,ind2) = X12e;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f_q_bc,q_bc,ind_bc1,sindq]=boundaryconditions_ConvectionVolume(globalnr_1h,globalnr_1v,XGLLGLL,YGLLGLL,velocity);

reorder_temp
break
Dp = T4*Dp;

ind_bc2 = [ 2*N2+(1:2:4*N) 2*N2+5*N+(1:2:2*N) ];
% Dp(ind_bc2,:) = [];

D_bc = Dp(:,ind_bc1);
Dp(:,ind_bc1)  = [];

XX = [X12 T3];
XX(ind_bc1,:) = [];
% XX(:,ind_bc2) = [];

A = Dp*XX;

A(ind_bc2,:) = [];
A(:,ind_bc2) = [];

D_bc(ind_bc2,:) = [];

f = -D_bc*f_q_bc;



% Remove boundary points

p = A\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% single_element_volume_postprocessen

etaGLL = xiGLL;
dxideta = diff(xiGLL)'*diff(etaGLL);

P = reshape(p(1:N2),N,N);

% P = P./dxideta;

for i=1:N
    for j=1:N
surf([xiGLL(i:i+1) ;  xiGLL(i:i+1)],[etaGLL(j) etaGLL(j) ; etaGLL(j+1) etaGLL(j+1)],P(i,j)/dxideta(i,j)*ones(2))
hold on
    end
end
xlabel('\xi')
ylabel('\eta')
view([0 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx = linspace(-1,1,200);
XX = xx'*ones(1,200);
YY = XX';

XiGLL  = xiGLL'*ones(1,N+1);
EtaGLL = ones(N+1,1)*etaGLL;

[hh,dhhdxx] = LagrangeVal(xx,N,1);
ee = EdgeVal(dhhdxx);
pp = ee'*P*ee;

figure
surf(XX,YY,pp)
% pcolor(XX,YY,pp)
shading interp
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 .6])
% set(gca,'clim',[-.1 1.1])
axis('square')