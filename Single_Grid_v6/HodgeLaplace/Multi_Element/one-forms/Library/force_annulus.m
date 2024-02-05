function F = force_annulus(R,Theta)

if length(Theta)~=5
    warning('wrong dimension Theta')
end
if length(R)~=2
    warning('wrong dimension R')
end

global N xi
global globalnr_1v globalnr_1h nr_1

eta = xi;

F = zeros(nr_1,1);
for nel=1:4
    
theta = (Theta(nel)+Theta(nel+1))/2 + (Theta(nel+1)-Theta(nel))/2*xi';
radius = (R(1)+R(2))/2 + (R(2)-R(1))/2*eta;

fv = zeros(N+1,N);
fh = zeros(N,N+1);

for i=1:N+1
    for j=1:N
        fv(i,j) = -1/2*cos(theta(i))^2*(radius(j+1)^2-radius(j)^2);
    end
end
for i=1:N
    for j=1:N+1
        fh(i,j) = 1/4*radius(j)^2*(cos(2*theta(i+1))-cos(2*theta(i)));
    end
end

F(globalnr_1v(:,nel)) = reshape(fv,N*(N+1),1);
F(globalnr_1h(:,nel)) = reshape(fh',N*(N+1),1);

end