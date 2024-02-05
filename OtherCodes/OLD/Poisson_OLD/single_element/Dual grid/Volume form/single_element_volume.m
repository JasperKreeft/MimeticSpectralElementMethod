clear all
% close all
% clc

path(path,'/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element/library_PoissonSingleElement')
path(path,'/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element/Dual grid/Library_DualGrid')

global N problem gridtype c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instellingen

method = 2;

NrCellRange = 12;

problem = 'sine' ;
gridtype = 'sinecurve' ;

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


switch gridtype
    case 'standard'
        c=0.0;
    case 'sinecurve' % from Hyman & Shaskov
        c=0.0; % Transformation in volumetozero not yet included
end

buildgrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method == 1
    [A,B,Dp,Gd,F] = elementmatrix_method1(N,xiGLL,xiEG,wG,wGLL,XiGLLGLL,EtaGLLGLL);
    
    %%%
    % Boundary conditions
    % since boundary conditions are zero (phi_bc=0),
    Dp = Dp(1:N^2,:);
    Gd = -Dp';
    
    C=volumetozero(N,xiG);

%     tic
%     qBphi = [ -A B*Gd*C ; Dp zeros(N^2) ]\[ zeros(2*N*(N+1),1) ; F];
%     phi_in = qBphi(2*N*(N+1)+1:end);
%     toc

    tic
    E = linsolve(full(A),full(B*Gd*C));
    phi_in = (Dp*E)\F;
    toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif method == 2
    [A,B,Dp,Gd,F] = elementmatrix_method2(N,xiGLL,xiEG,wG,wGLL,XiGLLGLL,EtaGLLGLL);

    %%%
    % Boundary conditions
    % since boundary conditions are zero (phi_bc=0),
    Dp = Dp(1:N^2,:);
    Gd = -Dp';

    C=volumetozero(N,xiG);
    
    tic
    qphi = [ A Dp'*B*C ; Dp zeros(N^2) ]\[ zeros(2*N*(N+1),1) ; F];
    phi_in = qphi(2*N*(N+1)+1:end);
    toc

%     tic
%     opts.SYM    = true;
%     opts.POSDEF = true;
%     opts.TRANSA = true;
% 
%     E = linsolve(full(A),full(Gd*B'*C),opts);
%     phi_in = (Dp*E)\F;
%     toc

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

single_element_volume_postprocessen

% fout_volume !!! Doet het nog niet !!!

end


% if length(NrCellRange)>1; errorplot; end !!! Doet het nog niet !!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element/library_PoissonSingleElement')
rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element/Dual grid/Library_DualGrid')
