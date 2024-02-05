clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%

global N N2 numRows numColumns RC
global nodes_in_element edges_in_element cells_in_element

numColumns = 5;
numRows    = 5;
RC = numColumns*numRows;

N = 5;

N2 = N*N;

a = 1;
b = 1;
velocity = [ a ; b ];

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);
cells_in_element = N2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering of 0-cells                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[XGLLGLL,YGLLGLL,XGG,YGG,QGLLGLL,JGLLGLL,JGG] = gridgenerator2();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dpe = topology(N);
Dpe = Dpe(1:N2,:);
% Divergence operator

Dp = kron(speye(RC),Dpe);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global xiGLL
global h e

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
e    = EdgeVal(dhdxi );


X12  = spalloc(RC*2*N*(N+1),RC*N2,RC*2*N2*(N+1));

for r=1:numRows
    for c=1:numColumns
        rc = (c+(r-1)*numColumns);

        % Convective part
        X12e = ElementConvectionVolumeMatrix(velocity);

        % Assembly
            ind1 = (rc-1)*edges_in_element+(1:edges_in_element);
            ind2 = (rc-1)*cells_in_element+(1:cells_in_element);
        X12(ind1,ind2) = X12e;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = numRows:-1:1
    for c = numColumns:-1:1
        max_in_element       = (c+(r-1)*numColumns)*edges_in_element;
        max_in_left_element  = ((c-1)+(r-1)*numColumns)*edges_in_element;
        max_in_lower_element = (c+(r-2)*numColumns)*edges_in_element;

        if r>1
            for i=1:N
                ind1 = max_in_lower_element-(i-1)*(N+1);
                ind2 = max_in_element-N-(i-1)*(N+1);
                Dp(:,ind1) = Dp(:,ind1) + Dp(:,ind2);
                Dp(:,ind2) = [];
                X12(ind2,:) = [];
            end
        end
        if c>1
            for i=1:N
                ind1 = max_in_left_element-N*(N+1)-(N+1)*(i-1);
                ind2 = max_in_element-(N+1)*(N+i)+1;
                Dp(:,ind1) = Dp(:,ind1) + Dp(:,ind2);
                Dp(:,ind2) = [];
                X12(ind2,:) = [];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f_q_bc,q_bc,ind_bc,sindq]=boundaryconditions_ConvectionVolume(globalnr_1h,globalnr_1v,XGLLGLL,YGLLGLL,velocity);

for c=2:numColumns
    ind_bc(c*N+1:end) = ind_bc(c*N+1:end)-N;
end
for r=3:numRows
    ind_bc((numColumns+r-1)*N+1:end) = ind_bc((numColumns+r-1)*N+1:end)-(2*numColumns-1)*N;
end

D_bc = Dp(:,ind_bc);

Dp(:,ind_bc)  = [];
X12(ind_bc,:) = [];

f = -D_bc*f_q_bc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indq = size(X12,1);
indp = size(X12,2);

B = [ -eye(indq) X12
      Dp zeros(indp) ];

F = [ zeros(indq,1)
          f         ];

qp = B\F;

q = qp(1:indq);
p = qp(indq+(1:indp));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

multi_element_volume_postprocessen
