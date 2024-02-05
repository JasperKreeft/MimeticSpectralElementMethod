clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N N2 numRows numColumns RC RCe
global nodes_in_element edges_in_element cells_in_element

numColumns = 4;
numRows    = 2;

N = 4;
N2 = N*N;

global m
m=1;

Pe = 0;
V = 0;

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
global wG

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


% velocity  = V*(1-YGLLGLL(1,:).^2);
% velocityG = V*(1-YGG(1,:).^2);

velocity  = V*ones(size((1-YGLLGLL(1,:).^2)));
velocityG = V*ones(size((1-YGG(1,:).^2)));


for r=1:numRows
    for c=1:numColumns
        rc = (c+(r-1)*numColumns);

        % Diffusive part
        [S1e,S0e,W0] = ElementSupportOperatorMethod(QGLLGLL(:,3*(rc-1)+(1:3)),JGLLGLL(:,rc));
%         [H02e]       = ElementHodgeMatrix(JGG(1:N2,rc),'02');
        [H02e]       = ElementHodgeMatrix(ones(N2,1),'02');

        S0e(1:N2,1:N2) = S0e(1:N2,1:N2)*H02e;

        % Convective part
        X12e = ElementConvectionVolumeMatrix(velocity((r-1)*N+(1:N+1)));

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

% A = [ S1/Re-eye(size(Dp,2))       Dp'*S0H02+X12
%        Dp   zeros(RC*(N2+2*N)+(numRows+numColumns)*N) ]; % FOUT

dT=1;
Told = 10;
qright = velocityG.*Told.*diff(YGLLGLL(1,:));
while dT>1e-3
    A = [ Pe*S1        Dp'*S0H02 + Pe*S1*X12
           Dp   zeros(RC*(N2+2*N)+(numRows+numColumns)*N) ];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Forcing function                                                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ind1 = RC*edges_in_element; % number of edges
    ind2 = RC*(N2+2*N)+(numRows+numColumns)*N; % number of points
    f1 = zeros(ind1,1);
    f2 = zeros(ind2,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Boundary conditions                                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ind_phi_bc1 = [globalnr_0(1,:) globalnr_0(2*(N+1)+(2:2*(N+1)),1)']; ind_phi_bc1(ind_phi_bc1==0)=[];

    ind_phi_bc2 = setdiff([globalnr_0(end,:) globalnr_0(:,end)' globalnr_0(2:2*(N+1),1)' globalnr_0(4*(N+1)+(2:(numColumns-4)*(N+1)),1)'],0);





    phi_bc	= [ ones(1,2*N)	2*ones(1,2*N) ]';
    qphi_bc	= zeros(N*2*(numColumns+numRows),1);
    qphi_bc([(numColumns-2)*N+(1:N) (2*numColumns-2)*N+(1:N)]) = 0*qright';

    ind_phi_bc1 = ind_phi_bc1+RCe;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_qbc = zeros(size(f2));
    f_pbc = zeros(size(f1));
    f_qbc(ind_phi_bc2,1) = qphi_bc(1:length(ind_phi_bc2));
    if sum(BC(2,:))>0
        f_pbc(1:RCe,1) = -A(1:RCe,ind_phi_bc1)*phi_bc;
    end

    f = [ f1+f_pbc
          f2+f_qbc ];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solving the system                                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A(:,ind_phi_bc1) = [];
    A(ind_phi_bc1,:) = [];
    f(ind_phi_bc1)   = [];

    qphi = A\f;

    q_in = qphi(1:RCe);
    phi_in = qphi(RCe+1:end);
    
    ind = [globalnr_0(end,2:N+1)-3*N globalnr_0(end,N+3:2*N+2)-4*N];
    T = phi_in(ind);
%     T = 1/2*T+1/2*Told; 
    uT = velocityG.*T';
    qright = velocityG.*T'.*diff(YGLLGLL(1,:));
    
    dT = sum(abs(T-Told));
    Told = T
break
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% postprocessen_q

postprocessen_me_volume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%