
n = N*Hconv+1;

E_ex = zeros(n*n,1);
for i=0:n-1
    for j=0:n-1
        ij = i+1+j*n;
        
        E_ex(ij,1) = i^2+j^2;
        
    end
end

E_ex = sort(E_ex);



plot(E_ex,'.r')
hold on