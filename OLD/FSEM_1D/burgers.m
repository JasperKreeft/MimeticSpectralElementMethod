% One-D Burgers

clear
clf
clc

Ne = 40;         % Nr of elements
Nn = Ne+1;      % Nr of nodes
L  = 1;         % length domain

% equidistant space
x = linspace(0,L,Ne+1);
h = diff(x);

bc_R = 0;
bc_L = 0;

% Connectivity matrix
c = zeros(Ne,2);
c(:,1) = (1:Ne)';
c(:,2) = (1:Ne)'+1;

M = zeros(Nn);
N = zeros(Nn);

for k=1:Ne
    
    Mass = [ h(k)/3  h(k)/6
             h(k)/6  h(k)/3];
    
    M(k:k+1,k:k+1) = M(k:k+1,k:k+1)+Mass;
    
end

% Periodic boundary condition
M([Nn 1],[Nn 1]) = M([Nn 1],[Nn 1])+Mass;
% N([Nn 1],[Nn 1]) = N([Nn 1],[Nn 1])+Adv;


% Solving the system
for i=1:Nn
    if x(i)<=0.5
        U(i,1) = 1/4-1/4*cos(4*pi*x(i));
    else
        U(i,1) = 0;
    end
end
plot(x,U)
t = 0;
dt = 1e-2;
Tend = 3;

while t<Tend
    N(1,[Nn 1 2]) = U(Nn)*[-1 1 0]/6+U(1)*[-1 0 1]/3+U(2)*[0 -1 1]/6;
    for i=2:Nn-1
        for j=i-1:i+1
            N(i,j) = U(i-1)*((i==j)-(i-1==j))/6+U(i)*((i+1==j)-(i-1==j))/3+U(i+1)*((i+1==j)-(i==j))/6;
        end
    end
    N(Nn,[1 Nn-1 Nn]) = U(Nn-1)*[-1 1 0]/6+U(Nn)*[-1 0 1]/3+U(1)*[0 -1 1]/6;

%     A = inv(M)*(M-dt*N); % Forward Euler
    A = inv(M+1/2*dt*N)*(M-1/2*dt*N); % Crank-Nicolson    
    t = t+dt;
    Unew = A*U;
    U = Unew;
    plot(x,U,'-o','linewidth',2)
    axis([0 1 -.1 .6])
    pause(0.01)
end