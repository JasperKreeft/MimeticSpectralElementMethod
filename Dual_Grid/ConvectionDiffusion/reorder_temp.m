T3e = spalloc(2*N*(N+1),4*N,4*N);
for i=1:N
    T3e((i-1)*(N+1)+1,2*i-1) = 1;
    T3e(i*(N+1),2*i) = 1;
    
    T3e((i-1+N)*(N+1)+1,2*(i+N)-1) = 1;
    T3e((i+N)*(N+1),2*(i+N)) = 1;
end
T3 = kron(speye(RC),T3e);

c=2;
for i=N:-1:1
    T3(:,(c-2)*4*N+2*i) = T3(:,(c-2)*4*N+2*i) + T3(:,(c-1)*4*N+2*i-1);
    T3(:,(c-1)*4*N+2*i-1) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = numColumns*(N2+3*N)+N;
T4 = spalloc(ind,ind,ind);

for c=1:numColumns
    T4((c-1)*N2+(1:N2),(c-1)*(N2+3*N)+(c>1)*N+(1:N2)) = speye(N2);
    T4(numColumns*N2+(c-1)*3*N+(c>1)*N+(1:3*N+(c==1)*N),c*N2+(c-1)*3*N+(c>1)*N+(1:3*N+(c==1)*N)) = speye(3*N+(c==1)*N);
end