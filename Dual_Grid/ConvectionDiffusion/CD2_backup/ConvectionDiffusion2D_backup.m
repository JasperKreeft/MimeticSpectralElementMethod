clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N N2 numRows numColumns RC
global nodes_in_element edges_in_element cells_in_element

numColumns = 2;
numRows    = 2;
RC = numColumns*numRows;

N = 2;
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

%      le ri lo up
BC = [ 1  0  0  0      % q
       0  1  1  1 ];   % phi

% ind = 2*N*(numRows+numColumns);
% qphi_bc = zeros(2*ind,1);


% indices

% Boundary indices for matrix A (all k-cells numbered)
ind_leftbc  = BC(1,1)*globalnr_1h(1,:)   ; ind_leftbc(ind_leftbc==0)   = [];
ind_rightbc = BC(1,2)*globalnr_1h(end,:) ; ind_rightbc(ind_rightbc==0) = [];
ind_lowerbc = BC(1,3)*globalnr_1v(:,1)'  ; ind_lowerbc(ind_lowerbc==0) = [];
ind_upperbc = BC(1,4)*globalnr_1v(:,end)'; ind_upperbc(ind_upperbc==0) = [];
ind_q_bc = sort([ind_leftbc ind_rightbc ind_lowerbc ind_upperbc]);

ind_leftbc  = BC(2,1)*globalnr_0(1,:)   ; ind_leftbc(ind_leftbc==0)   = [];
ind_rightbc = BC(2,2)*globalnr_0(end,:) ; ind_rightbc(ind_rightbc==0) = [];
ind_lowerbc = BC(2,3)*globalnr_0(:,1)'  ; ind_lowerbc(ind_lowerbc==0) = [];
ind_upperbc = BC(2,4)*globalnr_0(:,end)'; ind_upperbc(ind_upperbc==0) = [];
ind_phi_bc = sort([ind_leftbc ind_rightbc ind_lowerbc ind_upperbc])+RC*edges_in_element;

ind_bc = [ind_q_bc ind_phi_bc];

% Boundary indices for qphi_bc (only k-cells numbered with boundary condition)
ind_bc_lower = [N+1:N*numColumns N*(numColumns+1)+(1:N)];
ind_bc_left  = [(1:N) reshape((ones(numRows-1,1)*(1:N)+(0:2*N:2*N*(numRows-2))'*ones(1,N)+N*(numColumns+2))',1,[])];
ind_bc_right = [N*numColumns+(1:N) reshape((ones(numRows-2,1)*(1:N)+(N:2*N:2*N*(numRows-2))'*ones(1,N)+N*(numColumns+2))',1,[]) 2*N*(numRows+numColumns-1)+(1:N)];
ind_bc_upper = [N*(2*numRows+numColumns-1)+(1:N*(numColumns-1)) N*(2*numRows+2*numColumns-1)+(1:N)];
indq   = [BC(1,1)*ind_bc_left BC(1,2)*ind_bc_right BC(1,3)*ind_bc_lower BC(1,4)*ind_bc_upper]; indq(indq == 0)   = [];
indphi = [BC(2,1)*ind_bc_left BC(2,2)*ind_bc_right BC(2,3)*ind_bc_lower BC(2,4)*ind_bc_upper]; indphi(indphi==0) = [];

%%%

if BC(1,1)==true
    q_leftbc = (cos(pi*XGLLGLL(1,1:end-1)).*(cos(pi*YGLLGLL(1,1:end-1))-cos(pi*YGLLGLL(1,2:end))))';
else
    q_leftbc = [];
end
if BC(1,2)==true
    q_rightbc = 0*XGLLGLL(1,1:end-1)';
else
    q_rightbc = [];
end
if BC(1,3)==true
    q_lowerbc = 0*YGLLGLL(1:end-1,1);
else
    q_lowerbc = [];
end
if BC(1,4)==true
    q_upperbc = 0*YGLLGLL(1:end-1,end);
else
    q_upperbc = [];
end
q_bc = [q_leftbc ; q_rightbc ; q_lowerbc ; q_upperbc];
qphi_bc(indq,1) = q_bc;

%%%

if BC(2,1)==true
    phi_leftbc =0*(1+YGG(1,:)')/2;% 0.5*ones(size(YGG(1,:)')); %
else
    phi_leftbc = [];
end
if BC(2,2)==true
    phi_rightbc = 0*(1+YGG(end,:)')/2;%0.5*ones(size(YGG(end,:)')); %
else
    phi_rightbc = [];
end
if BC(2,3)==true
    phi_lowerbc = 0*0*XGG(:,1); % 1-XGG(:,1).^2;%0.5*ones(size(XGG(:,1)));%
else
    phi_lowerbc = [];
end
if BC(2,4)==true
    phi_upperbc = 0*ones(size(XGG(:,end)));%*XGG(:,end);
else
    phi_upperbc = [];
end
phi_bc = [phi_leftbc ; phi_rightbc ; phi_lowerbc ; phi_upperbc];
qphi_bc(indphi,1) = phi_bc;

%%%%

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

ind = RC*edges_in_element-size(indq,2);
q_in = qphi(1:ind);
phi_in = qphi(ind+1:end);
break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen_q

postprocessen_me

errorL2

errorL2_interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%