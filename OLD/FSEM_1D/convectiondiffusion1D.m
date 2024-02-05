% One-D convection-diffusion

clear
clf
clc

Ne = 10;         % Nr of elements
Nn = Ne+1;      % Nr of nodes
L  = 1;         % length domain

% equidistant space
% x = linspace(0,L,Ne+1);
% one-sided refinement
ratio = 1/100;
x = elm_line1(0,L,Ne,ratio);
h = diff(x);

Pe = 100;

bc_R = 0;
bc_L = 1;

% Connectivity matrix
c = zeros(Ne,2);
c(:,1) = (1:Ne)';
c(:,2) = (1:Ne)'+1;

F = zeros(Nn,1);
D = zeros(Nn);
N = zeros(Nn);

for k=1:Ne
    
    Adv  = [-1/2  1/2
            -1/2  1/2];
        
    Diff = [-1/h(k)  1/h(k)
             1/h(k) -1/h(k)];
    
    D(k:k+1,k:k+1) = D(k:k+1,k:k+1)+Diff;
    N(k:k+1,k:k+1) = N(k:k+1,k:k+1)+Adv;
    
end

K = -1/Pe*D+N;

% Left & right Dirichlet boundary condition
F(2) = F(2)-bc_L*K(2,1);
F(Nn-1) = F(Nn-1)-bc_R*K(Nn-1,Nn);

% Solving the system
% U = K(2:Nn-1,2:Nn-1)\F(2:Nn-1);
K_l = zeros(Nn,1); K_u = zeros(Nn,1);
K_d = diag(K); K_l(1:Nn-1) = diag(K,-1); K_u(2:Nn) = diag(K,1);
U = tridiag_band(K_l(2:Nn-1),K_d(2:Nn-1),K_u(2:Nn-1),F(2:Nn-1));

% Adding Dirichlet boundary condtions
U = [bc_L;U;bc_R];

figure(1)
plot(x,U,'-o','linewidth',2)
grid on
xlim([0 L])

hold on
xx = linspace(0,L,1000);
yy = bc_L+(bc_R-bc_L)*(exp(xx*Pe)-1)/(exp(Pe)-1);
plot(xx,yy,'r')