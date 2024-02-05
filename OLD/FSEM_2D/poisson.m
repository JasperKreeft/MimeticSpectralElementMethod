% Single element Galerkin SEM method for solving the poisson equation

clear all
close all
clc

N = 6;
m = 1;

[xGLL,wGLL] = GLLnodes(N);
XGLL = xGLL'*ones(1,N+1); YGLL = XGLL';

[h,dhdx] = LagrangeVal(xGLL,N,1);

A = zeros((N+1)^2);
for k=1:N+1
    for l=1:N+1
        kl = k+(l-1)*(N+1);
        for i=1:N+1
            for j=1:N+1
                ij = i+(j-1)*(N+1);
                A(kl,ij) = 2*sum(wGLL.*dhdx(i,:).*dhdx(k,:))+2*sum(wGLL.*dhdx(j,:).*dhdx(l,:));
            end
        end
    end
end

f = -2*m^2*pi^2*sin(m*pi*XGLL).*sin(m*pi*YGLL);

B = zeros((N+1)^2,1);
for p=1:N+1
    for q=1:N+1
        pq = p+(q-1)*(N+1);
        B(pq,1) = f(p,q)*wGLL(p)*wGLL(q);
    end
end

inner =[];
for j=2:N
    for i=2:N
        inner = [inner i+(j-1)*(N+1)];
    end
end

A = A(inner,inner);

B = B(inner,1);

u = A\B;