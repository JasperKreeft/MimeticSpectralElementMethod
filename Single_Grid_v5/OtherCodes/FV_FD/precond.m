
Nrange = [5 10 20 40];

k=0;
for N=Nrange
k=k+1;

load(['Afv_N' num2str(N) '.mat']);

Afv{k} = (A+A')/2;
clear A

load(['Amsemd_N' num2str(N) '.mat']);

Amsem{k} = (A+A')/2;
clear A

end

for i=1:k
    
    B = Afv{i}\Amsem{i};
    C = eig((B+B)/2);
    subplot(1,k,i)
    plot(abs(real(C)),imag(C),'.k')%
    axis equal
    axis([0 3 -1 1])
    title(['N = ' num2str(Nrange(i))])
    xlabel('Re')
    ylabel('Im')
    
end