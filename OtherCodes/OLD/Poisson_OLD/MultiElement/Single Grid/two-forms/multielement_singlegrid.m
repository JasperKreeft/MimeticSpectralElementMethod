clear all
close all
clc

global N m numRows numColumns
global xi w
global cc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings & loop for h-convergence                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
er = 0;

HconvRange = 2;

for Hconv = HconvRange

numRows    = 3;
numColumns = 3;
RC = numRows*numColumns;
NrCellRange = 2;
m = 1;
cc = 0.0;

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

localnr = reshape(1:N2,N,N);

globalnr_2 = zeros(numColumns*N,numRows*N);
for j=1:numRows
    for i=1:numColumns
        globalnr_2((i-1)*N+(1:N),(j-1)*N+(1:N)) = ((i-1)+(j-1)*numColumns)*N2+localnr;
    end
end

% disp('  '); disp('globalnr_2 = '); disp('  '); disp(num2str(flipud(globalnr_2')));

globalnr_1h = zeros(numColumns*N,numRows*N+1);
globalnr_1v = zeros(numColumns*N+1,numRows*N);
ind = 0;
for r=1:numRows
    for c=1:numColumns
        globalnr_1v((c-1)*N+(c>1)+(1:N+(c==1)),(r-1)*N+(1:N)) = ind + reshape(1:N*(N+(c==1)),N+(c==1),N) ;
        ind = ind + N2 + (c==1)*N;
        globalnr_1h((c-1)*N+(1:N),(r-1)*N+(r>1)+(1:N+(r==1))) = ind + reshape(1:N*(N+(r==1)),N+(r==1),N)';
        ind = ind + N*(N+(r==1));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);    eta = xi;  % Gauss-Lobotto-Legendre

[h,e] = MimeticpolyVal(xi,N,1);

Xi = xi'*ones(1,N+1); Eta = Xi';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numRows+1)   ;   % This might be more advanced !!!!!

% Global mesh
X = zeros(N*numColumns+1,N*numRows+1);
Y = zeros(N*numColumns+1,N*numRows+1);

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

D = zeros(globalnr_2(end,end),globalnr_1h(end,end));
ind = globalnr_1h(end,end);
M1 = zeros(ind);
ind1 = N2*numColumns*numRows;
ind2 = N2^2*numColumns*numRows;
M2 = spalloc(ind1,ind1,ind2);
ind = N2*numColumns*numRows;
f = zeros(ind,1);

for r=1:numRows
    for c=1:numColumns

        % Everything is in GLL-GLL nodes
        Xib  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xi )'*ones(1,N+1);
        Etab = ones(N+1,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*eta );

        % gridtype = 'sinecurve' % This might be more advanced !!!!!
        Xe = Xib + cc*sin(pi*Xib).*sin(pi*Etab);
        Ye = Etab+ cc*sin(pi*Xib).*sin(pi*Etab);

        if c<numColumns && r<numRows
            X((c-1)*N+(1:N),(r-1)*N+(1:N)) = Xe(1:N,1:N);
            Y((c-1)*N+(1:N),(r-1)*N+(1:N)) = Ye(1:N,1:N);
        elseif c<numColumns && r==numRows
            X((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = Xe(1:N,1:N+1);
            Y((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = Ye(1:N,1:N+1);
        elseif c==numColumns && r<numRows
            X((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = Xe(1:N+1,1:N);
            Y((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = Ye(1:N+1,1:N);
        elseif c==numColumns && r==numRows
            X((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = Xe;
            Y((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = Ye;
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gridtype = 'sinecurve'
dXdXib  = 1+pi*cc*cos(pi*Xib).*sin(pi*Etab);
dXdEtab = pi*cc*sin(pi*Xib).*cos(pi*Etab);
dYdXib  = pi*cc*cos(pi*Xib).*sin(pi*Etab);
dYdEtab = 1+pi*cc*sin(pi*Xib).*cos(pi*Etab);

dXibdXi   = (Xib(N+1,1)-Xib(1,1))/2*ones(N+1);
dXibdEta  = zeros(N+1);
dEtabdXi  = zeros(N+1);
dEtabdEta = (Etab(1,N+1)-Etab(1,1))/2*ones(N+1);

dXdXi  = dXdXib.*dXibdXi +dXdEtab.*dEtabdXi ;
dXdEta = dXdXib.*dXibdEta+dXdEtab.*dEtabdEta;
dYdXi  = dYdXib.*dXibdXi +dYdEtab.*dEtabdXi ;
dYdEta = dYdXib.*dXibdEta+dYdEtab.*dEtabdEta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = dXdXi.*dYdEta-dXdEta.*dYdXi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qinv11 = kron(reshape(( dXdXi./J),1,(N+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEta./J),1,(N+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEta./J),1,(N+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXi./J),1,(N+1)^2),[1 0])';

Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M1e = innerproduct_oneforms(e,J,Qinv);
ind1 = [ reshape(globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N)),N*(N+1),1)
         reshape(globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1))',N*(N+1),1) ];
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

M2e = innerproduct_twoforms(e,J);
ind2 = globalnr_2((c-1)*N+(1:N),(r-1)*N+(1:N));
M2(ind2,ind2) = M2e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Divergence operator

D(ind2,ind1) = D(ind2,ind1) + Dpe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        F = force(Xib(:,1),Etab(1,:));
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

[pu_bc,ind_u_bc,ind_p_bc] = boundaryconditions(globalnr_1h,globalnr_1v,globalnr_2,X,Y);

% ind_p_bc = sort(ind_p_bc);
% Matrix(ind_p_bc,:) = [];
% Matrix(:,ind_p_bc) = [];


f_qbc = zeros(size(f));
for i=1:4 % [ 1=Le 2=Ri 3=Lo 4=Up ]
    ind = ( ((i>1)+(i>2))*numRows + ((i>3)+(i>4))*numColumns )*N + ...
          (1:( ((i==1)+(i==2))*numRows + ((i==3)+(i==4))*numColumns )*N);
    f_qbc(ind_u_bc(ind)) = f_qbc(ind_u_bc(ind)) + pu_bc(ind);
end

ind = 2*N*(N+1) + ((numColumns-1)+(numRows-1))*(2*N2+N) + (numColumns-1)*(numRows-1)*2*N2 + ...
      - length(ind_p_bc)*0 ;

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
    semilogy(NrCellRange,errorL2)
    hold on
    semilogy(NrCellRange,errorL2_interp,'--r')
end
if length(HconvRange)>1
    loglog(2./(HconvRange),errorL2)
    hold on
    loglog(2./(HconvRange),errorL2_interp,'--xr')
end
errorL2

errorL2_interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%