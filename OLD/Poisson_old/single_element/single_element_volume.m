clear all
% close all
% clc

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
        c=0.0;
end

buildgrid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if method == 1
    [A,B,Dp,Gd,F] = elementmatrix_method1(N,xiGLL,xiEG,wG,wGLL,XiGLLGLL,EtaGLLGLL);
    tic
    C = linsolve(full(A),full(B*Gd));
    phi_in = (Dp*C)\F;
    toc
%     clear phi_in
%     tic
%     qBphi = [ -A B*Gd ; Dp zeros(N^2) ]\[ zeros(2*N*(N+1),1) ; F];
%     phi_in = qBphi(2*N*(N+1)+1:end);
%     toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif method == 2
    [A,B,Dp,Gd,F] = elementmatrix_method2(N,xiGLL,xiEG,wG,wGLL,XiGLLGLL,EtaGLLGLL);

    [C]=volumetozero(N,xiG);
    
%     tic
%     opts.SYM    = true;
%     opts.POSDEF = true;
%     opts.TRANSA = true;
% 
%     C = linsolve(full(A),full(Gd*B'),opts);
%     phi_in = (Dp*C)\F;
%     toc
%     clear phi_in
    tic
    qphi = [ A Dp'*B*C ; Dp zeros(N^2) ]\[ zeros(2*N*(N+1),1) ; F];
    phi = qphi(2*N*(N+1)+1:end);
    toc
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% single_element_postprocessen_volume

% fout_volume

end


if length(NrCellRange)>1; errorplot; end

