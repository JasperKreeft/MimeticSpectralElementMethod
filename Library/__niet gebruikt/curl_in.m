function Cp = curl_in(N)

nr_cells = N*N;

Cp_u = spdiags([ones(nr_cells,1) -ones(nr_cells,1)],[0 N],nr_cells,(N+1)*N);

unit = [-1 zeros(1,N-1) +1];

Unit = zeros(N,N*N+1);
for i=1:N
    ind = ((i-1)*N+1):(i*N+1);
    Unit(i,ind) = unit;
end
Cp_v = zeros(nr_cells,N*(N+1));
for j=1:N
    ind1 = (j-1)*N+(1:N);
    ind2 = j:N*N+j;
    Cp_v(ind1,ind2) = Unit;
end

Cp = sparse([Cp_u Cp_v]);