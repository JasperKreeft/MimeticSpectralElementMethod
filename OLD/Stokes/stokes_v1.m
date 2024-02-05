clear all
close all
clc

Re = 1;

nu = 1/Re;

m=1;

Z = 10;

N = Z

[xp,wGLL] = GLLnodes(N);
[xd,wG] = Gnodes(N);
xd_ex = [-1 xd 1];

yp = xp; yd = xd; yd_ex = xd_ex;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xp,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xp,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xd_ex,N,1);
[hw_w,dhwdxiw] = LagrangeVal(xd_ex,N,3);

e    = EdgeVal(dhdxi);
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);
ew_w = EdgeVal(dhwdxiw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


W1_11 = zeros(N*(N+1)); W2_11 = zeros(N*(N+1));
for l=1:N
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                W1_11(kl,ij) = ew(i,k)*e_w(l,j+1)*wGLL(k)*wG(j);
            end
        end
    end
end
for p=1:N+1
    for l=1:N
        for j=1:N
            pl = p+(l-1)*(N+1);
            pj = p+(j-1)*(N+1);
            W2_11(pl,pj) = wGLL(p)*sum(wG.*e_w(j,2:N+1).*e_w(l,2:N+1));
        end
    end
end

W1_22 = zeros(N*(N+1)); W2_22 = zeros(N*(N+1));
for l=1:N+1
    for k=1:N
        kl = k+(l-1)*N;
        for j=1:N+1
            for i=1:N
                ij  = i+(j-1)*N;
                W1_22(kl,ij) = e_w(k,i+1)*ew(j,l)*wG(i)*wGLL(l);
            end
        end
    end
end
for q=1:N+1
    for k=1:N
        for i=1:N
            qk = k+(q-1)*N;
            qi = i+(q-1)*N;
            W2_22(qk,qi) = wGLL(q)*sum(wG.*e_w(i,2:N+1).*e_w(k,2:N+1));
        end
    end
end

W1 = [W1_11 zeros(N*(N+1)); zeros(N*(N+1)) W1_22];
W2 = [W2_11 zeros(N*(N+1)); zeros(N*(N+1)) W2_22];

H1 = inv(W2)*W1;
% H1 = [inv(W2_11)*W1_11 zeros(N*(N+1)); zeros(N*(N+1)) inv(W2_22)*W1_22];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W3 = zeros((N+1)^2);
for l=1:N+1
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        for j=1:N+1
            for i=1:N+1
                ij = i+(j-1)*(N+1);
                W3(kl,ij) = wGLL(k)*wGLL(l)*ew(i,k)*ew(j,l);
            end
        end
    end
end
W4 = zeros((N+1)^2);
for q=1:N+1
    for p=1:N+1
        pq = p+(q-1)*(N+1);
        W4(pq,pq) = wGLL(p)*wGLL(q);
    end
end

H2 = inv(W4)*W3;
% H2 = ............ snellere oplossing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W5_12 = zeros(N*(N+1));
for q=1:N
    for k=1:N+1
        kq = k+(q-1)*(N+1);
        for j=1:N
            for p=1:N+1
                pj = p+(j-1)*(N+1);
                W5_12(kq,pj) = wGLL(p)*wG(q)*ew(k,p)*e_w(j,q);
            end
        end
    end
end
W5_21 = zeros(N*(N+1));
for l=1:N+1
    for p=1:N
        pl = p+(l-1)*N;
        for q=1:N+1
            for i=1:N
                iq = i+(q-1)*N;
                W5_21(pl,iq) = wG(p)*wGLL(q)*e_w(i,p)*ew(l,q);
            end
        end
    end
end
W6_11 = zeros(N*(N+1));
for q=1:N
    for k=1:N+1
        kq = k+(q-1)*(N+1);
        for i=1:N+1
            iq = i+(q-1)*(N+1);
            W6_11(kq,iq) = wG(q)*sum(wGLL.*ew(i,:).*ew(k,:)); % sum(wG.*ew_w(i,2:N+1).*ew_w(k,2:N+1));
        end
    end
end
W6_22 = zeros(N*(N+1));
for l=1:N+1
    for p=1:N
        pl = p+(l-1)*N;
        for j=1:N+1
            pj = p+(j-1)*N;
            W6_22(pl,pj) = wG(p)*sum(wGLL.*ew(j,:).*ew(l,:));
        end
    end
end

W5 = [zeros(N*(N+1)) W5_12; W5_21 zeros(N*(N+1))];
% W5 = [W5_12 zeros(N*(N+1)); zeros(N*(N+1)) W5_21];
% W5 = [W5_21 zeros(N*(N+1)); zeros(N*(N+1)) W5_12];
W6 = [W6_11 zeros(N*(N+1)); zeros(N*(N+1)) W6_22];

H3 = inv(W6)*W5;
% H3 = .................

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dp,Gd,NGp,Cd,Cp,Dd,Gp] = topology(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions for single element flow!!!
unit2 = [-1 0; zeros(N-1,2); 0 1];

Cd_bc_u = [eye(N+1) zeros(N+1); zeros((N-1)*(N+1),2*(N+1)); zeros(N+1) -eye(N+1)];

Cd_bc_v = kron(eye(N+1),unit2);

Cd_bc = [Cd_bc_u Cd_bc_v];


uv_bc = zeros(4*(N+1),1);
uv_bc(N+2:2*(N+1),1) = ones(N+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nuHdH = nu*H3*NGp*H2;

A = [   Dp*H1  zeros(N^2) ;
      nuHdH*Cd    Gd     ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = [    zeros(N*N,1)    ;
     -nuHdH*Cd_bc*uv_bc ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = A\F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = B(1:N*(N+1));
v = B(N*(N+1)+1:2*N*(N+1));
p = B(2*N*(N+1)+1:end);

U = reshape(u,N+1,N);
U = [uv_bc(1:N+1) U uv_bc(N+2:2*(N+1))];
V = reshape(v,N,N+1);
V = [uv_bc(2*(N+1)+1:2:4*(N+1))'; V; uv_bc(2*(N+1)+2:2:4*(N+1))'];
P = reshape(p,N,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

post_stokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%