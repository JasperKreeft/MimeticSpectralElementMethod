clear all
close all
clc

path(path,'/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element/library_PoissonSingleElement')

global N problem gridtype c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instellingen

NrCellRange = 12%2:2:20;

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
        c=0.2;
end

buildgrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A,B,Dp,Gd,F] = elementmatrix_singlegrid(N,xiGLL,xiEG,wG,wGLL,XiGLLGLL,EtaGLLGLL);
B = B;

%%%
% Boundary conditions
% since boundary conditions are zero (phi_bc=0),
Dp = Dp(1:N^2,:);
Gd = -Dp';

% tic
% qBphi = [ A Dp' ; Dp zeros(N^2) ]\[zeros(2*N*(N+1),1) ; F];
% q = qBphi(1:2*N*(N+1));
% phi_in = B\qBphi(2*N*(N+1)+1:end);
% toc

tic
qphi = [ A Dp'*B ; B*Dp zeros(N^2) ]\[zeros(2*N*(N+1),1) ; B*F];
q = qphi(1:2*N*(N+1));
phi_in = qphi(2*N*(N+1)+1:end);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

single_element_postprocessen_singlegrid

% fout_singlegrid

end


if length(NrCellRange)>1; errorplot; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path,'/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element/library_PoissonSingleElement')
