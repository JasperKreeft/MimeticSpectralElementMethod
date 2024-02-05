function [ E10, E21 ] = Incidence(N)

nn = (N+1)*(N+1);   % Number of nodes in the mesh
ne = 2*N*(N+1);     % Number of edges in the mesh
nv = N*N;           % Number of volumes in the mesh

row = zeros(4*N*(N+1),1);
col = zeros(4*N*(N+1),1);
val = zeros(4*N*(N+1),1);

counter = 1;
for j=1:N
    for i=1:N+1
        
        edge = i + (j-1)*(N+1);
        
        node = i + (j-1)*(N+1);
        row(counter) = edge;
        col(counter) = node;
        val(counter) = -1;
        
        counter = counter + 1;
        
        node = i + j*(N+1);
        row(counter) = edge;
        col(counter) = node;
        val(counter) = 1;
        
        counter = counter + 1;
        
    end
end

for j=1:N+1
    for i=1:N
        
        edge = N*(N+1) + i + (j-1)*N;
        
        node = i + (j-1)*(N+1);
        row(counter) = edge;
        col(counter) = node;
        val(counter) = 1;
        
        counter = counter + 1;
        
        node = i + 1 + (j-1)*(N+1);
        row(counter) = edge;
        col(counter) = node;
        val(counter) = -1;
        
        counter = counter + 1;
        
    end
end

E10 = sparse(row,col,val,ne,nn,counter-1);

row = zeros(4*N*N,1);
col = zeros(4*N*N,1);
val = zeros(4*N*N,1);


counter = 1;
for j=1:N
    for i=1:N
        volm = i + (j-1)*N;
        edge = i + (j-1)*(N+1);
        row(counter) = volm;
        col(counter) = edge;
        val(counter) = -1;
        
        counter = counter + 1;
        
        edge = i + 1 + (j-1)*(N+1);
        row(counter) = volm;
        col(counter) = edge;
        val(counter) = 1;
        
        counter = counter + 1;
        
        volm = i + (j-1)*N;
        edge = N*(N+1) + i + (j-1)*N;
        row(counter) = volm;
        col(counter) = edge;
        val(counter) = -1;
        
        counter = counter + 1;
        
        edge = N*(N+1) + i + j*N;
        row(counter) = volm;
        col(counter) = edge;
        val(counter) = 1;
        
        counter = counter + 1;
    end
end

E21 = sparse(row,col,val,nv,ne,counter-1);

% E21*E10