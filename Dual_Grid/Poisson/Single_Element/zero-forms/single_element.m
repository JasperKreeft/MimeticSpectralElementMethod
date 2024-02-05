clear all
close all
clc

% addpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/library_PoissonSingleElement')
% addpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/Dual grid/Library_DualGrid')
% path(path,'/media/My Passport/MSEM/MSEM_codes/Dual_Grid/Poisson/Single_Element/Library_SingleElement')
% path(path,'/media/My Passport/MSEM/MSEM_codes/Dual_Grid/Poisson/Single_Element/zero-forms/Library_ZeroForms')
path(path,'G:\MSEM\MSEM_codes\Dual_Grid\Poisson\Single_Element\zero-forms\Library_ZeroForms')
path(path,'G:\MSEM\MSEM_codes\Dual_Grid\Poisson\Single_Element\Library_SingleElement')
global N cc problem
global xi w wG
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instellingen

method = 2;
plot_figures  = 0;
error_figures = 1;

NrCellRange = 40;%2:2:20;

problem = 'sine' ;

cc = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problem
    case 'sine'
        global m                                                 %#ok<TLEV>
        m = 2; % number of waves
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

numbering('single')

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

F = reshape(force(xiGLL',xiGLL),nr_2,1);

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

    tic
    L = [ M1       Dp'*W20
          W20'*Dp spalloc(nr_2,nr_2,0) ];
    b = [zeros(nr_1,1) ; W20'*F];
    qphi   = L\b;
    q      = qphi(1:nr_1);
    phi_in = qphi(nr_1+1:nr_1+nr_2);
    toc
keyboard
%     % snelste voor uniforme grids
%     tic
%     L = [ M1       Dp'
%           Dp spalloc(nr_2,nr_2,0) ];
%     b = [zeros(nr_1,1) ; F];
%     qBphi = L\b;
%     q = qBphi(1:nr_1);
%     phi_in = W20\qBphi(nr_1+1:nr_1+nr_2);
%     toc

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
% PHI = reshape(phi_in,N,N);
% phi = hG'*PHI*hG;
PHI = [zeros(N+2,1) [zeros(1,N); reshape(phi_in,N,N); zeros(1,N)] zeros(N+2,1)];
phi = hEG(2:N+1,:)'*PHI(2:N+1,2:N+1)*hEG(2:N+1,:)+... % inner part
      hEG(1,:)'*PHI(1,2:N+1)*hGEG+...    % left side
      hEG(N+2,:)'*PHI(N+2,2:N+1)*hGEG+...      % right side
      hGEG'*PHI(2:N+1,1)*hEG(1,:)+...    % lower side
      hGEG'*PHI(2:N+1,N+2)*hEG(N+2,:);         % upper side

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

% filename = ['Pconv_H1_c' num2str(cc) '.mat'];
% save(filename,'N','NrCellRange','errorL2','errorL2_interp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/library_PoissonSingleElement')
% rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/Dual grid/Library_DualGrid')
% rmpath('/media/My Passport/MSEM/MSEM_codes/Dual_Grid/Poisson/Single_Element/Library_SingleElement')
% rmpath('/media/My Passport/MSEM/MSEM_codes/Dual_Grid/Poisson/Single_Element/zero-forms/Library_ZeroForms')