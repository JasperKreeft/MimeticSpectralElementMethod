function C = curl_out_3D(N)

%%

C11 = zeros(N*N*(N+1),N*(N+1)*(N+1));
for i=1:N

unit1 = kron(speye(N),kron(speye(N+1),[ zeros(1,i-1) -1 zeros(1,N-i)]));
unit2 = kron(speye(N),kron(speye(N+1),[ zeros(1,i-1) +1 zeros(1,N-i)]));

unit = [ unit1 spalloc(N*(N+1),N*(N+1),0) ] + [ spalloc(N*(N+1),N*(N+1),0) unit2 ];
     
C11((1:N*(N+1))+(i-1)*N*(N+1),:) = unit;
end
C11 = sparse(C11);

%%

C12 = zeros(N*N*(N+1),N*(N+1)*(N+1));

for i=1:N

unit1 = kron(speye(N+1),[ zeros(1,i-1) -1 zeros(1,N-i)]);
unit2 = kron(speye(N+1),[ zeros(1,i-1) +1 zeros(1,N-i)]);

for j=1:N
    C12((1:N+1)+(i-1)*(N+1)+(j-1)*N*(N+1),(1:N*(N+1))+(j-1)*N*(N+1)) = unit1;
    C12((1:N+1)+(i-1)*(N+1)+(j-1)*N*(N+1),(1:N*(N+1))+j*N*(N+1))     = unit2;
end

end
C12 = sparse(C12);

%%

C13 = spalloc(N*N*(N+1),N*(N+1)*(N+1),0);

%%

C21 = zeros(N*N*(N+1),N*(N+1)*(N+1));

for j=1:N
for i=1:N+1
    
    el = zeros(N+1,N); el(i,j) = 1;
    unit1 = kron(speye(N),+el);
    unit2 = kron(speye(N),-el);
    
    C21((1:N*(N+1))+(j-1)*N*(N+1),(1:N*(N+1))+(i-1)*N*(N+1)) = [ unit1 zeros(N*(N+1),N) ] + [ zeros(N*(N+1),N) unit2 ];
    
end
end
C21 = sparse(C21);

%%

C22 = C13;

%%

C23 = C12;

%%

C31 = C13;

%%

C32 = C21;

%%

C33 = zeros(N*N*(N+1),N*(N+1)*(N+1));

for i=1:N+1
    
    el = zeros(N+1,1); el(i,1) = 1;
    
    unit1 = kron(speye(N),kron(speye(N),+el));
    unit2 = -unit1;

    C33(:,(1:N*(N+1))+(i-1)*N*(N+1)) = [ unit1 zeros(N*N*(N+1),N) ] + [ zeros(N*N*(N+1),N) unit2 ];
    
end

C33 = sparse(C33);

%%

C = [ C11 C12 C13
      C21 C22 C23
      C31 C32 C33 ];
