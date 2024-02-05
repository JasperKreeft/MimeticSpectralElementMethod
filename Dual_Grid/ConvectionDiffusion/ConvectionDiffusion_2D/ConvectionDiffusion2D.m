clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N N2 numRows numColumns RC
global nodes_in_element edges_in_element cells_in_element

numColumns = 4;
numRows    = 4;
RC = numColumns*numRows;

N = 4;
N2 = N*N;

global m
m=1;

Re = 1;
velocity = [ 0 ; 0 ];

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

    ind1 = edges_in_element*RC;
    ind2 = edges_in_element^2*RC;
S1 = spalloc(ind1,ind1,ind2);
    ind1 = (nodes_in_element-N)*RC;
    ind2 = (nodes_in_element-N)^2*RC;
S0 = spalloc(ind1,ind1,ind2);
    ind1 = nodes_in_element*RC;
    ind2 = edges_in_element*RC;
    ind3 = (2*N)*(N+1)*(N+2)*RC;
X01  = spalloc(ind1,ind2,ind3);

for r=1:numRows
    for c=1:numColumns
        rc = (c+(r-1)*numColumns);

        % Diffusive part
        [S1e,S0e,W0] = ElementSupportOperatorMethod(QGLLGLL(:,3*(rc-1)+(1:3)),JGLLGLL(:,rc));

        % Convective part
        X01e = ElementConvectionMatrix(velocity,T1);

        % Assembly
            ind1 = (rc-1)*nodes_in_element+(1:nodes_in_element);
            ind2 = (rc-1)*edges_in_element+(1:edges_in_element);
        S1(ind2,ind2)  = S1e;
        S0(ind1,ind1)  = S0e;
        X01(ind1,ind2) = X01e;
    end
end

W0 = kron(speye(RC),W0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X01 = T2*X01;

for r = numRows:-1:1
    for c = numColumns:-1:1
        max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
        % max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
        % max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;

        if r>1
            ind = max_in_element-(2*(1:N)-1);
            S0(:,ind) = [];
            S0(ind,:) = [];
            W0(:,ind) = [];
            W0(ind,:) = [];
        end
        if c>1
            ind = max_in_element-2*N-(2*(1:N)-1);
            S0(:,ind) = [];
            S0(ind,:) = [];
            W0(:,ind) = [];
            W0(ind,:) = [];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct matrix                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = [ S1/Re    Dp'*S0
      S0'*Dp W0*X01*Gd ];

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
        i = (c-1)*N+(1:N);
        j = (r-1)*N+(1:N);
        ind = globalnr_2(i,j)+...
              +(r-1)*(2*numColumns+1)*N+...
              +(r==1)*(c>1)*(c-1)*N+(r>1)*numColumns*N+...
              +(c-1)*2*N+(c>1)*N;
        f2(ind) = F;

    end
end

f = [     f1
      S0'*f2 ];

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
phi_in = qphi(ind+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen_q
break
postprocessen_me

errorL2

errorL2_interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%