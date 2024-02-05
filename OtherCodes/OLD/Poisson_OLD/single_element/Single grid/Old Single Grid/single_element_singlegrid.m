clear all
close all
clc
% tic
global N problem gridtype c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instellingen

NrCellRange = 2:2:20;

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

[A,B,Dp,Gd,F] = elementmatrix_method2_singlegrid(N,xiGLL,xiEG,wG,wGLL,XiGLLGLL,EtaGLLGLL);


%%%
% Boundary conditions
% since boundary conditions are zero (phi_bc=0),
Dp = Dp(1:N^2,:);
Gd = -Dp';

%%%
%     tic
%     Bphi_in = -(Dp*inv(A)*Dp')\F;
%     phi_in = B\Bphi_in;
%     toc

%     % snelste voor uniforme grids
    tic
    qBphi = [ A Dp' ; Dp zeros(N^2) ]\[zeros(2*N*(N+1),1) ; F];
    q = qBphi(1:2*N*(N+1));
    phi_in = B\qBphi(2*N*(N+1)+1:end);
    toc

%     % snelste voor non-uniforme grids
%     tic
%     opts.SYM    = true;
%     if c<=0.3
%     opts.POSDEF = true;
%     opts.TRANSA = true;
%     end
% 
%     C = linsolve(full(A),full(Gd*B));%,opts);
%     phi_in = (Dp*C)\F;
%     q = C*phi_in;
%     toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% single_element_postprocessen_singlegrid

fout_singlegrid

end


if length(NrCellRange)>1; errorplot; end

% toc