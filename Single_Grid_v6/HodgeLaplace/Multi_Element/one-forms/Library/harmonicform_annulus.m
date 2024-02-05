function H = harmonicform_annulus(R,Theta)

if length(Theta)~=5
    warning('wrong dimension Theta')
end
if length(R)~=2
    warning('wrong dimension R')
end

global N xi
global globalnr_1v globalnr_1h nr_1

eta = xi;

H = zeros(nr_1,1);
for nel=1:4
    
theta = (Theta(nel)+Theta(nel+1))/2 + (Theta(nel+1)-Theta(nel))/2*xi';
radius = (R(1)+R(2))/2 + (R(2)-R(1))/2*eta;

hv = zeros(N+1,N);
hh = zeros(N,N+1);

for i=1:N+1
    for j=1:N
        hv(i,j) = log(radius(j+1))-log(radius(j));
        %(radius(j+1)*log(radius(j+1))-radius(j+1))-(radius(j)*log(radius(j))-radius(j));
    end
end
for i=1:N
    for j=1:N+1
        hh(i,j) = 0;
    end
end

H(globalnr_1v(:,nel)) = reshape(hv,N*(N+1),1);
H(globalnr_1h(:,nel)) = reshape(hh',N*(N+1),1);

end