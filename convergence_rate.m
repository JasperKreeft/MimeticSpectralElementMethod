clear all
close all
clc

kk = 1000;
xx = linspace(-1,1,kk);

kk0 = 500;
X0 = kk0/kk*2-1;

P = 8;

exact = zeros(1,kk);
for p=0:P-1
exact(1:kk0)    = exact(1:kk0)    + 1/factorial(p)*xx(1:kk0).^p;
end

for p=0:P
exact(kk0+1:kk) = exact(kk0+1:kk) + 1/factorial(p)*(1+(p~=P)*1/factorial(P-p)*(-X0)^(P-p))*xx(kk0+1:kk).^p;
end


plot(xx,exact)

Nrange = [1:6 8:2:32];

for N=Nrange

xgll = GLLnodes(N);
values = zeros(1,N+1);
for i=1:N+1
if xgll(i)<X0
for p=0:P-1
values(i) = values(i) + 1/factorial(p)*xgll(i)^p;
end
else
for p=0:P
values(i) = values(i) + 1/factorial(p)*(1+(p~=P)*1/factorial(P-p)*(-X0)^(P-p))*xgll(i)^p;
end
end
end

[hh,ee] = MimeticpolyVal(xx,N,1);

Interp = values*hh;

% hold on
% plot(xx,Interp,'r')


L2error(N) = sqrt(1/kk*sum((exact-Interp).^2));

end

figure
subplot(1,2,1)
semilogy(Nrange,L2error(Nrange))
subplot(1,2,2)
loglog(Nrange,L2error(Nrange))

%%

H = 10;
N = 1;

xigll = GLLnodes(N);

xixi = linspace(-1,1,kk/H+H);
[hh,ee] = MimeticpolyVal(xixi,N,1);

for h=1:H
    
    xgll = (xigll+1)/H+(2*(h-1)/H-1);
    values = zeros(1,N+1);
    for i=1:N+1
    if xgll(i)<X0
    for p=0:P-1
    values(i) = values(i) + 1/factorial(p)*xgll(i)^p;
    end
    else
    for p=0:P
    values(i) = values(i) + 1/factorial(p)*(1+(p~=P)*1/factorial(P-p)*(-X0)^(P-p))*xgll(i)^p;
    end
    end
    end
    
    kkh =  (h-1)*kk/H+(0:kk/H);
    
    interp(kkh(2:end)) = values*hh;
    
end






