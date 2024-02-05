clear all
close all
clc

addpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/library_PoissonSingleElement')
addpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/Dual grid/Library_DualGrid')

global N c problem
global xi w wG
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instellingen

method = 2;
plot_figures  = 1;
error_figures = 0;

NrCellRange = 12%:4:30;

problem = 'sine' ;

c = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problem
    case 'sine'
        global m                                                 %#ok<TLEV>
        m = 1; % number of waves
    case 'cosine'

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorL2 = zeros(size(NrCellRange));
ConditionNumber = zeros(size(NrCellRange));

for N=NrCellRange

disp(['N = ' num2str(N)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numbering

buildgrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h  ,e  ] = MimeticpolyVal(xiGLL,N,1);
[hw ,ew ] = MimeticpolyVal(xiGLL,N,3);  % only for method 1
[h_w,e_w] = MimeticpolyVal(xiEG ,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w = wGLL;
M1 = innerproduct_oneforms(e,JGLLGLL,Qinv);

if method == 1
    global wGL
    wGL = wGLL;
    W11 = wedgeproduct_1P1Dform(ew,e_w);
elseif method == 2
    W20 = wedgeproduct_2P0Dform(e_w);
end

Dp = topology(N);
% Boundary conditions
% since boundary conditions are zero (phi_bc=0)
Dp = Dp(1:N^2,:);

F = reshape(force(xiGLL,xiGLL),nr_2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method == 1
%         tic
%     L = [  M1   W11*Dp'
%            Dp zeros(N^2) ];
%     b = [ zeros(nr_1,1) ; F];
%        
%     qphi = L\b;
%     q = qphi(1:nr_1);
%     phi_in = qphi(nr_1+1:nr_1+nr_2);
%     toc

    % Snelste methode !!!!
    tic
    C = linsolve(full(M1),full(-W11*Dp'));
    phi_in = (Dp*C)\F;
    q = C*phi_in;
    toc
    
elseif method == 2
    %     tic
%     Bphi_in = -(Dp*inv(M1)*Dp')\F;
%     phi_in = W20\Bphi_in;
%     toc

%     tic
%     L = [ M1       Dp'*W20
%           W20'*Dp spalloc(nr_2,nr_2,0) ];
%     b = [zeros(nr_1,1) ; W20'*F];
%     qphi   = L\b;
%     q      = qphi(1:nr_1);
%     phi_in = qphi(nr_1+1:nr_1+nr_2);
%     toc

    % snelste voor uniforme grids
    tic
    L = [ M1       Dp'
          Dp spalloc(nr_2,nr_2,0) ];
    b = [zeros(nr_1,1) ; F];
    qBphi = L\b;
    q = qBphi(1:nr_1);
    phi_in = W20\qBphi(nr_1+1:nr_1+nr_2);
    toc

%     % snelste voor non-uniforme grids
%     tic
%     opts.SYM    = true;
%     if c<=0.3
%     opts.POSDEF = true;
%     opts.TRANSA = true;
%     end

%     C = linsolve(full(M1),full(Gd*W20));%,opts);
%     phi_in = (Dp*C)\F;
%     q = C*phi_in;
%     toc
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSEN                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid, basis-functions and weights for post-processen
postproces_grid

% Exact Solution
exact_solution


% Reconstruct zero-form
% Add boundary conditions to solution
% phi = reconstruct_zeroforms(PHI,hEG)
% PHI = [zeros(N+2,1) [zeros(1,N); reshape(phi_in,N,N); zeros(1,N)] zeros(N+2,1)];
% phi = hEG'*PHI*hEG;
PHI = reshape(phi_in,N,N);
phi = hG'*PHI*hG;


% Gradient
Gd = -Dp';
grad_phi = Gd*phi_in;
U = grad_phi(1:N*(N+1),1);
V = grad_phi(N*(N+1)+1:2*N*(N+1),1);

u = [zeros(N+1,1) reshape(U,N+1,N) zeros(N+1,1)];
v = [zeros(1,N+1); reshape(V,N+1,N)'; zeros(1,N+1)];

uxi      = eEG'*u*hEG;
ueta     = hEG'*v*eEG;
ux = (dYdEtap.*uxi-dYdXip.*ueta)./Jp;
uy = (dXdXip.*ueta-dXdEtap.*uxi)./Jp;

velo = sqrt(ux.^2+uy.^2);

% Reconstruction of Flux
[qx,qy,qq] = reconstruct_oneforms(q,hGL,eGL,Jp,dXdXip,dXdEtap,dYdXip,dYdEtap);

% Plotten
if plot_figures
    plotten
end

% error
if error_figures
    fout
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % for N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(NrCellRange)>1 && error_figures
    errorplot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/library_PoissonSingleElement')
rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/Dual grid/Library_DualGrid')
