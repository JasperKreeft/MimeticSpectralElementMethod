clear all
close all
clc

for N = 1:20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADJOINT

[xgl,wgl] = GLLnodes(N);

[h,e]=MimeticpolyVal(xgl,N,1);

D = full(topology1D(N));

M0 = diag(wgl);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(wgl.*e(i,:).*e(j,:));
    end
end

Dstar = inv(M0)*D'*M1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -star d star
[xg,wg] = Gnodes(N);
xeg = [-1 xg 1]; weg = [0 wg 0];

[hw,ew]   = MimeticpolyVal(xgl,N,3);
[h_w,e_w] = MimeticpolyVal(xg,N,1);

H1 = e_w';

G = -D';

H2 = ew';

Dstar2 = -H2*G*H1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L(N) = norm(sum(Dstar,2)-sum(Dstar2,2));
LL(N) = norm(sum(Dstar,1)-sum(Dstar2,1));

end