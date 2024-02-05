clear all
close all
clc

global N m numRows numColumns
global xi
global cc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings & loop for h-convergence                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er = 0;

HconvRange = 5;

for Hconv = HconvRange

numRows    = Hconv;
numColumns = Hconv;
RC = numRows*numColumns;
NrCellRange = 5%1:2:14;
m = 1;
cc = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for convergence plots                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(Hconv)
% 
% errorL2        = zeros(1,max(NrCellRange));
% errorL2_interp = zeros(1,max(NrCellRange));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop for p-convergence                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for N = NrCellRange

disp(['N = ' num2str(N)]);

N2 = N*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering unknowns                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering_singlegrid();

nr_1 = globalnr_1h(end,end);
nr_2 = globalnr_2(end,end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,Qinv,J] = gridgenerator_singlegrid();

eta = xi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the divergence matrix for the internal points.
Dpe = zeros(N2,2*N*(N+1));
for i = 1:N
    for j = 1:N
        Dpe(i + (j-1)*N,i+(j-1)*(N+1))   = -1 ;
        Dpe(i + (j-1)*N,i+1+(j-1)*(N+1)) =  1 ;
        Dpe(i + (j-1)*N,N*(N+1) + j+(i-1)*(N+1))     = -1 ;
        Dpe(i + (j-1)*N,N*(N+1) + j+1+(i-1)*(N+1))   =  1 ;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = zeros(nr_2,nr_1);
M1 = zeros(nr_1);
M2 = zeros(nr_2);
f = zeros(nr_2,1);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

ind1 = [ reshape(globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N)),N*(N+1),1)
         reshape(globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1))',N*(N+1),1) ];
ind2 = globalnr_2((c-1)*N+(1:N),(r-1)*N+(1:N));

% Divergence operator
D(ind2,ind1) = Dpe;
     
% one-forms
Qinve = spdiags(Qinv(:,3*(rc-1)+(1:3)),-1:1,2*(N+1)^2,2*(N+1)^2);
M1e = innerproduct_oneforms(e,J(:,rc),Qinve);
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct_twoforms(e,J(:,rc));
M2(ind2,ind2) = M2e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        F = force_twoform(r,c);

        ind = globalnr_2((c-1)*N+(1:N),(r-1)*N+(1:N));
        f(ind) = F;

    end
end

% pcolor(X,Y,zeros(size(X))); axis equal; axis([-1 1 -1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct matrix                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = N2*numColumns*numRows;

Matrix = full([ M1  D'*M2 ; M2'*D spalloc(ind,ind,0) ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary Conditions & Right-hand-side                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [pu_bc,ind_u_bc,ind_p_bc] = boundaryconditions(globalnr_1h,globalnr_1v,globalnr_2,X,Y);
% 
% f_qbc = zeros(size(f));
% for i=1:4 % [ 1=Le 2=Ri 3=Lo 4=Up ]
%     ind = ( ((i>1)+(i>2))*numRows + ((i>3)+(i>4))*numColumns )*N + ...
%           (1:( ((i==1)+(i==2))*numRows + ((i==3)+(i==4))*numColumns )*N);
%     f_qbc(ind_u_bc(ind)) = f_qbc(ind_u_bc(ind)) + pu_bc(ind);
% end
% 
% ind = 2*N*(N+1) + ((numColumns-1)+(numRows-1))*(2*N2+N) + (numColumns-1)*(numRows-1)*2*N2 + ...
%       - length(ind_p_bc)*0 ;

ind = 2*N*(N+1) + ((numColumns-1)+(numRows-1))*(2*N2+N) + (numColumns-1)*(numRows-1)*2*N2 ;

RHS = [  zeros(ind,1)
               M2'*f ]; %       M2'*(f+f_qbc) ];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pu = Matrix\RHS;

ind = numel(globalnr_1h) + numel(globalnr_1v);
% ind = numel(globalnr_1h(:,2:end-1)) + numel(globalnr_1v(2:end-1,:));

p = pu(1:ind,1);
u = pu(ind+(1:numRows*numColumns*N2),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen_me_singlegrid

end % for N
end % for H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for convergence plots                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% errorL2

if length(NrCellRange)>1
    figure
    semilogy(NrCellRange,errorL2)
    hold on
    semilogy(NrCellRange,errorL2_interp,'--r')
    grid on
end
if length(HconvRange)>1
    figure
    loglog(2./(HconvRange),errorL2)
    hold on
    loglog(2./(HconvRange),errorL2_interp,'--xr')
    grid on
end
errorL2

errorL2_interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%