clear all
close all
clc

N = 2;

x = linspace(-1,1,1000);

[xp,wp] = GLLnodes(N);
[xd,wd] = Gnodes(N);

[hpd,epd] = MimeticpolyVal(xd,N,1);
[hdp,edp] = MimeticpolyVal(xp(2:N),N,2);

D = topology1D(N);
D = D(:,2:N);
Dd = -D';

Dstar1 = -edp'*Dd*epd'

%%%%%%%%%%

[hpp,epp] = MimeticpolyVal(xp,N,1);

M0 = diag(wp(2:N));

for i=1:N
    for j=1:N
        M1(i,j) = sum(wp.*epp(i,:).*epp(j,:));
    end
end

Dstar2 = inv(M0)*D'*M1