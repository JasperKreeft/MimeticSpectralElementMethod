clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings & loop for h-convergence                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 0;
for Hconv = 2

numRows    = 2;%Hconv;
numColumns = 2;%Hconv;
NrCellRange = 3;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% standard element mesh and weights
[xiGLL,wGLL] = GLLnodes(N);  etaGLL = xiGLL;
[xiG,wG]     = Gnodes(N);    etaG   = xiG;
xiEG         = [-1. xiG 1.]; etaEG  = xiEG;

% Main mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]
x = linspace(-1,1,numColumns+1);
y = linspace(-1,1,numRows+1);


% Global mesh
XGLL = zeros(N*numColumns+1,N*numRows+1);
YGLL = zeros(N*numColumns+1,N*numRows+1);
XG   = zeros(N*numColumns  ,N*numRows  );
YG   = zeros(N*numColumns  ,N*numRows  );
for j=1:numRows
    for i=1:numColumns

        XeGLL = ( (x(i+1)+x(i))/2+(x(i+1)-x(i))/2*xiGLL )'*ones(1,N+1);
        YeGLL = ones(N+1,1)*( (y(j+1)+y(j))/2+(y(j+1)-y(j))/2*etaGLL );
        XeG   = ( (x(i+1)+x(i))/2+(x(i+1)-x(i))/2*xiG )'*ones(1,N);
        YeG   = ones(N,1)*( (y(j+1)+y(j))/2+(y(j+1)-y(j))/2*etaG );

        if i<numColumns && j<numRows
            XGLL((i-1)*N+(1:N),(j-1)*N+(1:N)) = XeGLL(1:N,1:N);
            YGLL((i-1)*N+(1:N),(j-1)*N+(1:N)) = YeGLL(1:N,1:N);
        elseif i<numColumns && j==numRows
            XGLL((i-1)*N+(1:N),(j-1)*N+(1:N+1)) = XeGLL(1:N,1:N+1);
            YGLL((i-1)*N+(1:N),(j-1)*N+(1:N+1)) = YeGLL(1:N,1:N+1);
        elseif i==numColumns && j<numRows
            XGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N)) = XeGLL(1:N+1,1:N);
            YGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N)) = YeGLL(1:N+1,1:N);
        elseif i==numColumns && j==numRows
            XGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N+1)) = XeGLL;
            YGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N+1)) = YeGLL;
        end
        XG((i-1)*N+(1:N),(j-1)*N+(1:N)) = XeG;
        YG((i-1)*N+(1:N),(j-1)*N+(1:N)) = YeG;

    end
end

xGLL = XGLL(:,1); yGLL = YGLL(1,:);

surf(XGLL,YGLL,zeros(size(XGLL))); view([0 0 1])
clear XG YG XGLL YGLL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dpe,Gde] = topology(N);

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

Dbc = zeros(N*numRows*(numColumns-1)+N*(numRows-1)*numColumns,numRows*numColumns*edges_in_element);
j = 0;
for r = 1:numRows
    for c = 1:numColumns
        max_in_element       = (c+(r-1)*numColumns)*edges_in_element;
        max_in_left_element  = ((c-1)+(r-1)*numColumns)*edges_in_element;
        max_in_upper_element = (c+r*numColumns)*edges_in_element;
        if c<numColumns
            for i=1:N
                j = j + 1;
                % horizontal fluxes
                Dbc( j , max_in_left_element+i*(N+1)  ) = -1;
                Dbc( j , max_in_element+1+(i-1)*(N+1) ) = +1;
            end
        end
        if r<numRows
            for i=1:N
                j = j + 1;
                % vertical fluxes
                Dbc( j , max_in_element-(N-i) ) = -1;
                Dbc( j , max_in_upper_element-N2-(N-i) ) = +1;
            end
        end
    end
end

D = [Dp ; Dbc ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
He = elementmatrix(N,xiGLL,xiEG,wG,wGLL,numRows,numColumns);

H = kron(speye(numRows*numColumns),He);

Q = H*G;        clear H;

A=D*Q;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% F = force_v2(m,XGLL(:,1),YGLL(1,:));
F = force_v2(m,xGLL,yGLL);

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

postprocessen_me

end % for N
end % for H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for convergence plots                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% errorL2
% 
% semilogy(NrCellRange,errorL2(NrCellRange))
% hold on
% semilogy(NrCellRange,errorL2_interp(NrCellRange),'--r')
% 
% loglog(2./(2:2:Hconv),errorL2)
% hold on
% loglog(2./(2:2:Hconv),errorL2_interp,'--xr')

errorL2

errorL2_interp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%