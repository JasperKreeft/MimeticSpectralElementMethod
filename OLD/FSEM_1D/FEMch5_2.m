clear
clf
clc

bc_R = 0;
bc_L = 0;

% Exact solution
xx = 0:0.01:1;
yy = (bc_L*exp(1)-bc_R+1)/(exp(1)-exp(-1))*exp(-xx)+(bc_R-bc_L*exp(-1)-1)/(exp(1)-exp(-1))*exp(xx)+xx;
plot(xx,yy,'g','linewidth',2)

tic
Ne = 2;         % Nr of elements
Nn = 2*Ne+1;      % Nr of nodes
L  = 1;         % length domain

% equidistant space
% x = linspace(0,L,Ne+1);
% one-sided refinement
ratio = 1/1;
x = elm_line1 (0,L,Ne,ratio);
h = diff(x);



elements = zeros(Ne,3);
elements(:,1) = (1:2:Nn-2)';
elements(:,2) = (2:2:Nn-1)';
elements(:,3) = (3:2:Nn)';

F = zeros(Nn,1);
K = zeros(Nn);
K_ll = zeros(Nn,1);
K_l  = zeros(Nn,1);
K_d  = zeros(Nn,1);
K_u  = zeros(Nn,1);
K_uu = zeros(Nn,1);

for k=1:Ne
    Kk11 = 2/15*h(k)+7/(3*h(k)); Kk12 = 1/15*h(k)-8/(3*h(k));  Kk13 = -1/30*h(k)+1/(3*h(k));
    Kk21 = Kk12;                 Kk22 = 8/15*h(k)+16/(3*h(k)); Kk23 = 1/15*h(k)-8/(3*h(k));
    Kk31 = Kk13;                 Kk32 = Kk23;                  Kk33 = 2/15*h(k)+7/(3*h(k));
%     
%     K_ll(2*k+1) = Kk31;
%     K_l(2*k)    = Kk21;
%     K_l(2*k+1)  = Kk32;
%     K_d(2*k-1)  = Kk11+K_d(2*k-1);
%     K_d(2*k)    = Kk22;
%     K_d(2*k+1)  = Kk33;
%     K_u(2*k-1)  = Kk12;
%     K_u(2*k)    = Kk23;
%     K_uu(2*k-1) = Kk13;

    K((2*k-1):(2*k+1),(2*k-1):(2*k+1)) = K((2*k-1):(2*k+1),(2*k-1):(2*k+1))+...
                                         [Kk11 Kk12 Kk13;
                                          Kk21 Kk22 Kk23;
                                          Kk31 Kk32 Kk33];

    Fk = h(k)/6*[x(k); 2*x(k)+2*x(k+1); x(k+1)];
    for i=1:3
        F(elements(k,i)) = F(elements(k,i))+Fk(i);
    end
end

F(2:3) = F(2:3)-bc_L*[Kk21;Kk31];

F(end-2:end-1) = F(end-2:end-1)-bc_R*[Kk13;Kk23];

A = K(2:end-1,2:end-1)\F(2:end-1);
% A = penta(Nn-2,K_d(2:end-1),K_u(2:end-1),K_uu(2:end-1),K_l(2:end-1),K_ll(2:end-1),F(2:end-1));
A = [bc_L;A;bc_R];

x_2 = (x(1:Ne)+x(2:Ne+1))/2;

X(1:2:Nn) = x;
X(2:2:Nn) = x_2;
toc

hold on
for k=1:Ne
YI = interp1(X((2*k-1):(2*k+1)),A((2*k-1):(2*k+1))',linspace(x(k),x(k+1),20),'spline');
plot(linspace(x(k),x(k+1),20),YI,'r','linewidth',2)
end


plot(X(2:2:Nn),A(2:2:Nn),'xk','linewidth',2)
plot(X(1:2:Nn),A(1:2:Nn),'sk','linewidth',2)
grid on