clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left boundary is periodic boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 20;

a = 2;

%%%%%%%%%%%%%%%%%%%%%%

[Dp,Gd,Cd,NGp,Cp,Dd,Gp] = topology(N);
clear Dp Gd Cd NGp Cp Dd
% per bc.
Gp(:,(N+1):(N+1):(N+1)^2) = Gp(:,(N+1):(N+1):(N+1)^2)+Gp(:,1:(N+1):(N+1)^2);
Gp(:,1:(N+1):(N+1)^2) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xiGLL = GLLnodes(N);

[h,dhdxi] = LagrangeVal(xiGLL,N,1);
e         = EdgeVal(dhdxi);

IT = kron(speye(N+1),e');

IT = [ a*IT zeros((N+1)^2,N*(N+1))
       zeros((N+1)^2,N*(N+1)) IT ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_of_nodes = (N+1)^2;
C1 = spalloc(nr_of_nodes,2*nr_of_nodes,2*nr_of_nodes);

C1(1:nr_of_nodes,1:nr_of_nodes) = speye(nr_of_nodes);
for j=1:N+1
    C1((j-1)*(N+1)+(1:N+1),nr_of_nodes+(j-1)+(1:N+1:nr_of_nodes)) = speye(N+1);
end
% per bc.
C1((N+1):(N+1):(N+1)^2,:) = C1((N+1):(N+1):(N+1)^2,:)+C1(1:(N+1):(N+1)^2,:);
C1(1:(N+1):(N+1)^2,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
A = C1*IT*Gp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_ic = zeros(1,N);
for i=1:N
    if xiGLL(i+1)<=0
        phi_ic(1,i) = 1/4-1/4*cos(2*pi*xiGLL(i+1));
    else
        phi_ic(1,i) = 0;
    end
end


B = A;
B(1:N,:) = [];
f=-B(:,1:N)*phi_ic';
B(:,1:N) = [];

Phi = B\f;

phi = reshape(Phi,N,N)';

phi = [phi_ic ; phi ];

phi = [phi(:,N) phi];

figure(1)
subplot(1,2,1)
pcolor(xiGLL,xiGLL,phi)
shading interp
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 .6])
axis('square')



xx = linspace(-1,1,201);

hh = LagrangeVal(xx,N,1);
pphi = hh'*phi*hh;

subplot(1,2,2)
pcolor(xx,xx,pphi)
shading interp
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 .6])
axis('square')

figure(2)
subplot(1,2,1)
plot([linspace(-1,0,100) 0 1],[1/4-1/4*cos(2*pi*linspace(-1,0,100)) 0 0],'--r')
hold on
plot(xx,pphi(end,:),xiGLL,phi(end,:),'ob','markerface','b')

for i=1:200
    figure(2)
    subplot(1,2,2)
    plot(xx,pphi(i,:))%,xiGLL,phi(i,:),'ob','markerface','b')
    pause(0.1)
end


figure(3)
plot([linspace(-1,0,100) 0 1],[1/4-1/4*cos(2*pi*linspace(-1,0,100)) 0 0],'--r')
hold on
plot(xx,pphi(101,:))
plot(xx,pphi(201,:),'g')