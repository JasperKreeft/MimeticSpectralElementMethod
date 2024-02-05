clear all
% close all
clc

global N m cc numRows numColumns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings & loop for h-convergence                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 0;

HconvRange = 2;

for Hconv = HconvRange

numRows    = Hconv;
numColumns = Hconv;
NrCellRange = 3;%:2:25;
m = 1;

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

globalnr = zeros(numColumns*N,numRows*N);
for j=1:numRows
    for i=1:numColumns
        globalnr((i-1)*N+(1:N),(j-1)*N+(1:N)) = ((i-1)+(j-1)*numColumns)*N2+localnr;
    end
end

disp('  '); disp('globalnr = '); disp('  '); disp(num2str(flipud(globalnr')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% grid option
cc = 0.0;

buildgrid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keyboard


[Dpe,Gde] = topology_new(N);


% Gradient operator

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);

Gd_bc_r = zeros(edges_in_element,4*N);
Gd_bc_l = zeros(edges_in_element,4*N);
Gd_bc_a = zeros(edges_in_element,4*N);
Gd_bc_b = zeros(edges_in_element,4*N);


for i=1:N
    Gd_bc_r(i*(N+1),N+2*i)                = +1;
    Gd_bc_l((i-1)*(N+1)+1,N-1+2*i)        = -1;
    Gd_bc_a(edges_in_element-(N-i),3*N+i) = +1;
    Gd_bc_b(N*(N+1)+i,i)                  = -1;
end

Gd_bc = Gd_bc_r+Gd_bc_l+Gd_bc_b+Gd_bc_a;

G = kron(speye(numRows*numColumns),[Gde Gd_bc]);

for r = numRows:-1:1
    for c = numColumns:-1:1
        max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
        max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
        max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;
        % horizontal elimination
        if c>1
            for i=1:N
                ind1 = max_in_left_element-N-2*(i-1);
                ind2 = max_in_element-N-(2*i-1);
                G( : , ind1 ) = G( : , ind1) + G( : , ind2 );
                G( : , ind2 ) = [];
            end
        end
        % vertical elimination
        if r>1
            for i=1:N
                ind1 = max_in_lower_element-(i-1);
                ind2 = max_in_element-3*N-(i-1);
                G( : , ind1 ) = G( :, ind1 )+ G( : , ind2 );
                G( : , ind2 ) = [];
            end
        end
    end
end

% Divergence operator

Dp = kron(speye(numRows*numColumns),Dpe);

for r = numRows:-1:1
    for c = numColumns:-1:1
        if r>1
            
            
        end
        
        if c>1
            
            
        end
        
        
    end
end












% Dbc = zeros(N*numRows*(numColumns-1)+N*(numRows-1)*numColumns,numRows*numColumns*edges_in_element);
% j = 0;
% for r = 1:numRows
%     for c = 1:numColumns
%         max_in_element       = (c+(r-1)*numColumns)*edges_in_element;
%         max_in_left_element  = ((c-1)+(r-1)*numColumns)*edges_in_element;
%         max_in_upper_element = (c+r*numColumns)*edges_in_element;
%         if c<numColumns
%             for i=1:N
%                 j = j + 1;
%                 % horizontal fluxes
%                 Dbc( j , max_in_left_element+i*(N+1)  ) = -1;
%                 Dbc( j , max_in_element+1+(i-1)*(N+1) ) = +1;
%             end
%         end
%         if r<numRows
%             for i=1:N
%                 j = j + 1;
%                 % vertical fluxes
%                 Dbc( j , max_in_element-(N-i) ) = -1;
%                 Dbc( j , max_in_upper_element-N2-(N-i) ) = +1;
%             end
%         end
%     end
% end
% Dbc = sparse(Dbc);
% 
% D = [Dp ; Dbc ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [Dpe,He,Gde] = elementmatrix_method1(N,xiGLL,xiEG,wG,wGLL,XibGLLGLL,EtabGLLGLL);
% H = kron(speye(numRows*numColumns),He);

% [Dpe,He,Gde] = elementmatrix_method2(N,xiGLL,xiEG,wG,wGLL,XibGLLGLL,EtabGLLGLL);
% Method 2 werkt nog niet, omdat He niet bestaat!!

H = elementmatrix_method1_curved(N,xiGLL,xiEG,wG,wGLL,xibLR,etabAB);

Q = H*G;        clear H;

A=D*Q;          clear Q;

if Hconv>1
% remove bc columns

max_in_element = (numRows*numColumns)*N2+(numRows*(numColumns-1)+numColumns*(numRows-1))*N+... % internal points
                 2*N*(numColumns+numRows); % boundary points
for r = numRows:-1:1
    for c = numColumns:-1:1
        if r == numRows
            for i=N:-1:1
%                 max_in_element-(N-i)
                A( : , max_in_element-(N-i) ) = [];
            end
        end
        if c == numColumns
            for i=N:-1:1
%                 max_in_element-N-(N-i)
                A( : , max_in_element-N-(N-i) ) = [];
            end
        end
        if c == 1
            for i=N:-1:1
%                 max_in_element-N-1-2*(N-i)
                A( : , max_in_element-N-1-2*(N-i) ) = [];
            end
        end
        if r == 1
            for i=N:-1:1
%                 max_in_element-(2+1*(c==1))*N-(N-i)
                A( : , max_in_element-(2+1*(c==1))*N-(N-i) ) = [];
            end
        end
        max_in_element = max_in_element-(N2+2*N)-N*((c==1)+(r==1));
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = force(XGLLGLL(:,1),YGLLGLL(1,:));
keyboard
f = zeros(N*numColumns*N*numRows,1);
for j=1:N*numRows
    for i=1:N*numColumns
        f(globalnr(i,j)) = F(i,j);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fbc = zeros((numRows*(numColumns-1)+numColumns*(numRows-1))*N,1);

f = [ f ; fbc];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = A\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HIERRR
% kiezen tussen G of EG of combinatie???
postprocessen_me_curved

end % for N
end % for H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for convergence plots                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% errorL2

if length(NrCellRange)>1
    semilogy(NrCellRange,errorL2(NrCellRange))
    hold on
    semilogy(NrCellRange,errorL2_interp(NrCellRange),'--r')
end
if length(HconvRange)>1
    loglog(2./(HconvRange),errorL2)
    hold on
    loglog(2./(HconvRange),errorL2_interp,'--xr')
end
errorL2

errorL2_interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%