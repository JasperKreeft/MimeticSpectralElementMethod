clear all; close all; clc;

N=2; N2 = N^2;

numRows    = 2;
numColumns = 2;

[Dpe,Gde] = topology_new(N);

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);
cells_in_element = N2;

% Divergence operator

Dp = kron(speye(numRows*numColumns),Dpe);

for r = numRows:-1:1
    for c = numColumns:-1:1
        max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
        max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
        max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;
        
        if r>1
            for i=1:N
                ind1 = max_in_lower_element-2*(i-1);
                ind2 = max_in_element-(2*i-1);
                Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
                Dp(ind2,:) = [];
            end            
        end
        if c>1
            for i=1:N
                ind1 = max_in_left_element-2*N-2*(i-1);
                ind2 = max_in_element-2*N-(2*i-1);
                Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
                Dp(ind2,:) = [];
            end
        end
    end
end

% Gradient operator
Gd = -Dp';