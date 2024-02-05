clear all
% clf%ose all
clc

N = 4;
H = 20;

Re = 100;

phiL = 1;
phiR = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact

nn = 100;

[xixi,ww] = GLLnodes(nn);
% xx = (0+1)/2+(1-0)/2*xixi;


Jac = 1/(2*H);
XiXi = zeros(1,nn*H+1);
for h=1:H
    ind = nn*(h-1)+(1:nn+1);
    XiXi(ind) = xixi+2*(h-1);
end
xx = (XiXi+1)*Jac;

phi_ex = (exp(Re*xx)-exp(Re))/(1-exp(Re));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

D = spdiags([-ones(H*N,1) ones(H*N,1)],[0 1],H*N,H*N+1);

E = zeros(N*H,N*H+1);
for h=1:H
    E((h-1)*N+(1:N),(h-1)*N+(1:N+1)) = e;
end
% Central
for h=1:H
    E(:,(h-1)*N+1) = E(:,(h-1)*N+1)/2;
end
% % Upwind
% for h=1:H-1
%     E(h*N+(1:N),h*N+1) = zeros(N,1);
% end
% % Downwind
% for h=1:H-1
%     E((h-1)*N+(1:N),h*N+1) = zeros(N,1);
% end
% % 1/2 Upwind, 1/2 central
% a = 5.175;      % a=2 is central     % Re = 100, N = 20, a = 5.175  % Re = 100, N = 20, a = 3.2;
% E(1:N,1)     = E(1:N,1)*(a-1)/a;
% E(N+(1:N),1) = E(N+(1:N),1)/a;
% for h=1:H-1
%     E((h-1)*N+(1:N),h*N+1) = E((h-1)*N+(1:N),h*N+1)*(a-1)/a;
%     E(h*N+(1:N),h*N+1) = E(h*N+(1:N),h*N+1)/a;
% end

Jac = 2;
M0 = zeros(1,N*H+1);
for h=1:H
    ind = N*(h-1)+(1:N+1);
    M0(ind) = M0(ind) + w*Jac;
end
M0 = diag(M0);

Jac = 2*(2*H);
M1e = zeros(N);
for i=1:N
    for j=1:N
        M1e(i,j) = sum(w.*e(i,:).*e(j,:)*Jac);
    end
end
M1 = kron(eye(H),M1e);


Conv = M0*E'*D;
Diff = -D'*M1*D;

Matrix = -Re*Conv + Diff;

% Boundary conditions

Matrix([1 N*H+1],:) = [];

Rhs = - phiL*Matrix(:,1) + phiR*Matrix(:,N*H+1);

Matrix(:,[1 N*H+1]) = [];

phi_in = Matrix\Rhs;

phi = [ phiL ; phi_in ; phiR ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jac = 1/(2*H);
Xi = zeros(1,N*H+1);
for h=1:H
    ind = N*(h-1)+(1:N+1);
    Xi(ind) = xi+2*(h-1);
end
x = (Xi+1)*Jac;


hh = LagrangeVal(xixi,N,1);

pphi = zeros(1,nn*H+1);
phi_interp = zeros(1,nn*H+1);
for h=1:H
    ind1 = nn*(h-1)+(1:nn+1);
    ind2 = N*(h-1)+(1:N+1);
    pphi(ind1) = phi(ind2)'*hh;
    phi_interp(ind1) = (exp(Re*x(ind2))-exp(Re))/(1-exp(Re))*hh;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure plot
plot(xx,phi_ex,'g')
hold on
plot(xx,pphi,'')
plot(x,phi,'o','markerface','b')
grid on
hold off
