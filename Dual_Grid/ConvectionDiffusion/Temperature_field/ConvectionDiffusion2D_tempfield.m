clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N N2 numRows numColumns RC RCe
global nodes_in_element edges_in_element cells_in_element

numColumns = 2;
numRows    = 2;

N = 4;
N2 = N*N;

global m
m=1;

Re = 1;
velocity = [ 0 ; 0 ];

global BC

%      le ri lo up
BC = [ 0  1  0  1      % q
       1  0  1  0 ];   % phi

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
S0 = spalloc(ind1,ind1,ind);

    ind = (2*N)*(N+1)*(N+2)*RC;
X01  = spalloc(ind1,ind2,ind);

    ind = RC*N2+((numRows-1)*numColumns+(numColumns-1)*numRows)*N;
H20W0 = spalloc(ind1,ind1,ind);

for r=1:numRows
    for c=1:numColumns
        rc = (c+(r-1)*numColumns);

        % Diffusive part
        [S1e,S0e,W0e] = ElementSupportOperatorMethod(QGLLGLL(:,3*(rc-1)+(1:3)),JGLLGLL(:,rc));
        [H20e]        = ElementHodgeMatrix(JGG(1:N2,rc),'20');

        % Convective part
        X01e = ElementConvectionMatrix(velocity,T1);

        % Assembly
            ind1 = (rc-1)*nodes_in_element+(1:nodes_in_element);
            ind2 = (rc-1)*edges_in_element+(1:edges_in_element);
            ind3 = (rc-1)*nodes_in_element+(1:cells_in_element);
        S1(ind2,ind2)  = S1e;
        S0(ind1,ind1)  = S0e;
        X01(ind1,ind2) = X01e;
        H20W0(ind1,ind1) = W0e;
        H20W0(ind3,ind3) = H20e*H20W0(ind3,ind3);
        
    end
end

W0 = kron(speye(RC),W0e);

X01 = T2*X01;
% break
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
            S0(:,ind2) = [];
            S0(ind2,:) = [];
            H20W0(:,ind2) = [];
            H20W0(ind2,:) = [];
            Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
            Dp(ind2,:) = [];
        end
        if c>1
            ind1 = max_in_left_element-2*N-2*((1:N)-1);
            ind2 = max_in_element-2*N-(2*(1:N)-1);
            S0(:,ind2) = [];
            S0(ind2,:) = [];
            H20W0(:,ind2) = [];
            H20W0(ind2,:) = [];
            Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
            Dp(ind2,:) = [];
        end
    end
end

% Gradient operator
Gd = -Dp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct matrix                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A = [ S1/Re    Dp'*S0
%       S0'*Dp W0*X01*Gd ];
%   
A = [  S1/Re   Dp'*S0
      S0'*Dp  H20W0*X01*Gd ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind1 = RCe; % number of edges
ind2 = RC*(N2+2*N)+(numRows+numColumns)*N; % number of points
f1 = zeros(ind1,1);
f2 = zeros(ind2,1);
for r=1:numRows
    for c=1:numColumns

%         F = force(r,c);
%         i = 1+(c-1)*(N+1)+(1:N);
%         j = 1+(r-1)*(N+1)+(1:N);
%         f2(globalnr_0(i,j)) = F;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[qphi_bc,q_bc,phi_bc,ind_q_bc1,ind_q_bc2,ind_phi_bc1,ind_phi_bc2,indq,indphi]=boundaryconditions(globalnr_0,globalnr_1h,globalnr_1v,XGLLGLL,YGLLGLL,XGG,YGG);
ind_phi_bc1 = ind_phi_bc1+RCe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_qbc = zeros(size(f2));
f_qbc(ind_phi_bc2,1) = qphi_bc(1:((BC(1,1)+BC(1,2))*numRows+(BC(1,3)+BC(1,4))*numColumns)*N);
f_pbc(1:RCe,1) = -A(1:RCe,ind_phi_bc1)*phi_bc;

f = [     (f1+f_pbc)
      S0'*(f2+f_qbc) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A(:,ind_phi_bc1) = [];
A(ind_phi_bc1,:) = [];
f(ind_phi_bc1)   = [];

qphi = A\f;

q_in = qphi(1:RCe);
phi_in = qphi(RCe+1:end);

if sum(BC(1,:))==4
    phi_in = phi_in - sum(phi_in)/length(phi_in);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen_q

postprocessen_me
view([0 0 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%