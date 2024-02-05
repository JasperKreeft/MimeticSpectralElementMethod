clear all
close all
clc

global N w e

N = 3;

[xi,w] = GLLnodes(N);
eta = xi;
zeta = xi;

[h,e] = MimeticpolyVal(xi,N,1);

G = grad_in_3D(N);

M0 = innerproduct_3D(0);
M1 = innerproduct_3D(1);

Matrix_full = G'*M1*G;
Matrix = Matrix_full;

F = zeros((N+1)^3,1);
for i=1:N+1
    for j=1:N+1
        for k=1:N+1
            ijk = i+(j-1)*(N+1)+(k-1)*(N+1)^2;
            F(ijk,1) = -3*pi^2*sin(pi*xi(i))*sin(pi*eta(j))*sin(pi*zeta(k));
        end
    end
end

RHS = M0*F;

%%

globalnr0 = zeros(N+1,N+1,N+1);
for k=1:N+1
    globalnr0(:,:,k) = (k-1)*(N+1)^2+reshape(1:(N+1)^2,N+1,N+1);
end

bottom = reshape(globalnr0(:,:,1),(N+1)^2,1);
top    = reshape(globalnr0(:,:,N+1),(N+1)^2,1);
front  = reshape(globalnr0(:,1,:),(N+1)^2,1);
back   = reshape(globalnr0(:,N+1,:),(N+1)^2,1);
left   = reshape(globalnr0(1,:,:),(N+1)^2,1);
right  = reshape(globalnr0(N+1,:,:),(N+1)^2,1);

bc_number = unique([ bottom
                     top
                     front
                     back
                     left
                     right ]);

%%

Matrix(bc_number,:) = [];
Matrix(:,bc_number) = [];
RHS(bc_number,:)    = [];

phi = Matrix\RHS;

%%

Phi_ex = zeros(N+1,N+1,N+1);
for i=1:N+1
    for j=1:N+1
        for k=1:N+1
            Phi_ex(i,j,k) = sin(pi*xi(i))*sin(pi*eta(j))*sin(pi*zeta(k));
        end
    end
end


%%

% Phi = zeros(N+1,N+1,N+1);
% for i=1:N+1
%     for j=1:N+1
%         for k=1:N+1
%             ijk = i+(j-1)*(N+1)+(k-1)*(N+1)^2;
%             Phi(i,j,k) = phi(ijk,1);
%         end
%     end
% end

Phi = zeros(N+1,N+1,N+1);
for i=2:N
    for j=2:N
        for k=2:N
            ijk = i-1+(j-2)*(N-1)+(k-2)*(N-1)^2;
            Phi(i,j,k) = phi(ijk,1);
        end
    end
end

Xi = xi'*ones(1,N+1);
Eta = ones(N+1,1)*eta;

%%

for k=1:N+1
     surf(Xi,Eta,Phi(:,:,k)+2*k)
    hold on
end
