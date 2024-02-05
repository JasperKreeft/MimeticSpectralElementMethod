clear all
close all
clc

N = 2;

u = 1;
v = 1;

[xgl,wgl] = GLLnodes(N); ygl = xgl;
[xg,wg]   = Gnodes(N);

[h,dhdx] = LagrangeVal(xgl,N,1);
e = EdgeVal(dhdx);

[ni nj] = size(e);
f = zeros(ni-1,nj);
for i=1:ni-1
    for j=1:nj
        f(i,j) = -sum(e(1:i,j));
    end
end




for k=1:N
    for l=1:N
        kl = k+(l-1)*N;
        for i=2:N
            for j=1:N
                ij = 2*j-1+(i-2)*N;
                U1(kl,ij)   = -sum(wgl.*f(i-1,:).*e(k,:))*sum(wgl.*e(j,:).*e(l,:));
                U1(kl,ij+1) = -U1(kl,ij);
            end
        end
    end
end

for k=1:N
    for l=1:N
        kl = k+(l-1)*N;
        for i=1:N
            for j=2:N
                ij = i+(j-2)*N;
                V1(kl,ij)   = -sum(wgl.*e(i,:).*e(k,:))*sum(wgl.*f(j-1,:).*e(l,:));
                V1(kl,ij+N) = -V1(kl,ij);
            end
        end
    end
end


for k=2:N
    for l=1:N
        kl = 2*l-1+(k-2)*N;
        for i=1:N
            for j=1:N
                ij = i+(j-1)*N;
                U2(kl,ij) = -sum(wgl.*e(i,:).*f(k-1,:))*sum(wgl.*e(j,:).*e(l,:));
                U2(kl+1,ij) = -U2(kl,ij);
            end
        end
    end
end

for k=1:N
    for l=2:N
        kl = k+(l-2)*N;
        for i=1:N
            for j=1:N
                ij = i+(j-1)*N;
                V2(kl,ij) = -sum(wgl.*e(i,:).*e(k,:))*sum(wgl.*e(j,:).*f(l-1,:));
                V2(kl+N,ij) = -V2(kl,ij);
            end
        end
    end
end

U = (U1-U2);
V = (V1-V2);

for k=1:N
    for l=1:N
        kl = k+(l-1)*N;
        for i=1:N
            for j=1:N
                ij = i+(j-1)*N;
                Ml(kl,ij) = e(i,1)*e(k,1)*sum(wgl.*e(j,:).*e(l,:));
                Mr(kl,ij) = e(i,N+1)*e(k,N+1)*sum(wgl.*e(j,:).*e(l,:));
                Mb(kl,ij) = sum(wgl.*e(i,:).*e(k,:))*e(j,1)*e(l,1);
                Ma(kl,ij) = sum(wgl.*e(i,:).*e(k,:))*e(j,N+1)*e(l,N+1);
            end
        end
    end
end

abcb = v*xgl+u;
abcr = -v-u*ygl;

for k=1:N
    for l=1:N
        kl = k+(l-1)*N;
        bcb(kl,1) = v*e(l,1)*sum(wgl.*abcb.*e(k,:));
        bcr(kl,1) = u*e(k,1)*sum(wgl.*abcr.*e(l,:));
    end
end

bc = 1/2*(bcb+bcr);



A = 1/2*(u*U+v*V) + 1/2*(u*Mr+v*Ma);

rhs = bc;

a = A\rhs;


nn = 100;
xx = linspace(-1,1,nn);
[hh,dhhdx] = LagrangeVal(xx,N,1);
ee = EdgeVal(dhhdx);

aa = ee'*reshape(a,N,N)*ee;
xx = xx'*ones(1,nn);
yy = xx';
pcolor(xx,yy,aa)
shading interp
axis equal
axis([-1 1 -1 1])
colorbar

