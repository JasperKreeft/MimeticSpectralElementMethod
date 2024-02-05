clear all
close all
clc

global N m numRows numColumns
global xi w
global cc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings & loop for h-convergence                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er = 0;

HconvRange = [2 4 8 16 32 ];

for Hconv = HconvRange

numRows     = Hconv;
numColumns  = Hconv;
NrCellRange = 2;
m = 1;
cc = 0.0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for convergence plots                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(Hconv)
% 
% errorL2        = zeros(1,max(NrCellRange));
% errorL2_interp = zeros(1,max(NrCellRange));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop for p-convergence                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for N = NrCellRange

disp(['N = ' num2str(N)]);

N2 = N*N;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering unknowns                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering_singlegrid();

nr_0 = globalnr_0(end,end);
nr_1 = globalnr_1h(end,end);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);    eta = xi;  % Gauss-Lobotto-Legendre

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,Qinv,J] = gridgenerator_singlegrid();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gpev = [ -eye(N*(N+1)) zeros(N*(N+1),N+1) ] + [ zeros(N*(N+1),N+1) eye(N*(N+1)) ];
Gpeu = zeros(N*(N+1),(N+1)^2);
for i=1:N
    Gpeu((i-1)*(N+1)+(1:N+1),:) = kron(eye(N+1),[ zeros(1,i-1) 1 -1 zeros(1,N-i) ]);
end

Gpe = [ Gpev ; Gpeu ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = zeros(nr_1,nr_0);
M0 = zeros(nr_0); %spalloc(nr_0,nr_0,nr_0);
M1 = zeros(nr_1);
f = zeros(nr_0,1);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

ind1 = reshape(globalnr_0((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)),1,[]);
ind2 = [ reshape(globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N)),N*(N+1),1)
         reshape(globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1))',N*(N+1),1) ];

% Gradient operator
G(ind2,ind1) = Gpe;
     
% zero-forms
M0e = innerproduct_zeroforms(J(:,rc));
M0(ind1,ind1) = M0(ind1,ind1) + M0e;

% one-forms
Qinve = spdiags(Qinv(:,3*(rc-1)+(1:3)),-1:1,2*(N+1)^2,2*(N+1)^2);
M1e = innerproduct_oneforms(e,J(:,rc),Qinve);
M1(ind2,ind2) = M1(ind2,ind2) + M1e;

M0 = sparse(M0);
M1 = sparse(M1);
G  = sparse(G);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Xe = X((c-1)*N+(1:N+1),(r-1)*N+(1:N+1));
        Ye = Y((c-1)*N+(1:N+1),(r-1)*N+(1:N+1));

        F = force_zeroform(Xe,Ye);
        ind = globalnr_0((c-1)*N+(1:N+1),(r-1)*N+(1:N+1));
        f(ind) = F;

    end
end

% figure
% pcolor(X,Y,zeros(size(X))); axis equal; axis([-1 1 -1 1]);
% 
% figure
% surf(X,Y,f(globalnr_0)); shading interp


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct matrix                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrix = -G'*M1*G;
RHS = M0*f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Element- and Domain boundary Conditions                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boundary_points = sort([ globalnr_0(:,1)' globalnr_0(:,end)' globalnr_0(1,2:end-1) globalnr_0(end,2:end-1) ]);

Matrix(boundary_points,:) = [];
Matrix(:,boundary_points) = [];

RHS(boundary_points) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = Matrix\RHS;

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
end
if length(HconvRange)>1
    figure
    loglog(2./(HconvRange),errorL2)
    hold on
    loglog(2./(HconvRange),errorL2_interp,'--xr')
end
errorL2

errorL2_interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%