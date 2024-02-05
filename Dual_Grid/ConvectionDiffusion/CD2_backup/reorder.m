function [T1,T2] = reorder()

global N N2 numRows numColumns 
global nodes_in_element edges_in_element cells_in_element

T1 = spalloc(nodes_in_element,2*(N2+2*N),2*(N2+2*N));

for j=1:N
    ind1 = (j-1)*N+(1:N);
    ind2 = (j-1)*(N+2)+1+(1:N);
    T1(ind1,ind2) = speye(N);
end

ind1 = N2+(1:2:2*N);
ind2 = 1:N+2:N*(N+2);
T1(ind1,ind2) = speye(N);

ind1 = N2+(2:2:2*N);
ind2 = (1:N)*(N+2);
T1(ind1,ind2) = speye(N);

for i=1:N
    ind1 = (1:N:N2)+(i-1);
    ind2 = N*(N+2)+1+(1:N)+(i-1)*(N+2);
    T1(ind1,ind2) = speye(N);
end

ind1 = N2+2*N+(1:2:2*N);
ind2 = N*(N+2)+(1:N+2:N*(N+2));
T1(ind1,ind2) = speye(N);

ind1 = N2+2*N+1+(1:2:2*N);
ind2 = N*(N+2)+(1:N)*(N+2);
T1(ind1,ind2) = speye(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T2 = speye(numRows*numColumns*(N2+4*N));
for r=numRows:-1:1
    for c=numColumns:-1:1

        max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
        max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
        max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;

        if r>1
            ind1 = max_in_left_element+N2+2*N+(1:2:2*N);
            ind2 = max_in_lower_element-(2*N-2:-2:0);
            T2(ind2,ind1) = speye(N);
            T2(ind1,:) = [];
        end

        if c>1
            ind1 = max_in_left_element+N2+(1:2:2*N);
            ind2 = max_in_left_element-4*N+(2:2:2*N);
            T2(ind2,ind1) = speye(N);
            T2(ind1,:) = [];
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%