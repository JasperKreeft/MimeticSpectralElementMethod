clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%

N = 20;
M = N;

a = .5;

%%%%%%%%%%%%%%%%%%%%%%

[Dp,Gd,Cd,NGp,Cp,Dd,Gp] = topology(N);
clear Dp Gd Cd NGp Cp Dd
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
A = C1*IT*Gp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


phi_ic = zeros(1,N+1);
for i=1:N+1
    if xiGLL(i)<=0
        phi_ic(1,i) = 1/4-1/4*cos(2*pi*xiGLL(i));
    else
        phi_ic(1,i) = 0;
    end
    
%     if xiGLL(i)>-0.4
%         phi_ic(1,i) = 1;
%     end
end

B = A;
B([1:N+1 (1:N)*(N+1)+1],:) = [];
f=-B(:,2:N+1)*phi_ic(2:N+1)';
B(:,[1:N+1 (1:N)*(N+1)+1]) = [];

Phi = B\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = reshape(Phi,N,N)';

phi = [phi_ic(2:N+1) ; phi ];

phi = [zeros(N+1,1) phi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,2,1)
pcolor(xiGLL,xiGLL,phi)
shading interp
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 .6])
% set(gca,'clim',[-.1 1.1])
axis('square')

xx = linspace(-1,1,200);

hh = LagrangeVal(xx,N,1);
pphi = hh'*phi*hh;

subplot(1,2,2)
pcolor(xx,xx,pphi)
shading interp
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 .6])
% set(gca,'clim',[-.1 1.1])
axis('square')

figure
subplot(1,2,1)
plot([-1 0 linspace(0,1,100)],[0 0 1/4-1/4*cos(2*pi*linspace(0,1,100))],'--r')
hold on
plot(xx,pphi(end,:),xiGLL,phi(end,:),'ob','markerface','b')


% pause(1)
% for i=1:200
%     figure(2)
%     subplot(1,2,2)
%     plot(xx,pphi(i,:))%,xiGLL,phi(i,:),'ob','markerface','b')
% %     axis([-1 1 0.45 0.55])
% %     grid on
%     pause(0.05)
% end

figure
surf(xx,xx,pphi,'FaceColor','interp',...
	'EdgeColor','none',...
	'FaceLighting','phong')
shading interp
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 0.6])
% set(gca,'clim',[-.1 1.1])
axis tight
view(-50,30)
camlight left




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Teruglopen

% IT = [ -a*IT_GLLGLL_s IT_GLLGLL_t ];
% 
% A = IT*G;
% 
% phi_ic = phi(end,:);
% 
% B = A;
% B([1:N+1 (1:N)*(N+1)+1],:) = [];
% f=(-B(:,2:N+1)*phi_ic(2:N+1));
% B(:,[1:N+1 (1:N)*(N+1)+1]) = [];
% 
% Phi = B\f;
% 
% phi = reshape(Phi,N,N)';
% 
% phi = [phi_ic(2:N+1)' ; phi ];
% 
% phi = [zeros(N+1,1) phi];
% 
% subplot(1,2,1)
% pcolor(xiGLL,xiGLL,phi)
% shading interp
% xlabel('x')
% ylabel('t')
% colorbar
% set(gca,'clim',[-.1 .6])
% axis('square')
% 
% 
% 
% xx = linspace(-1,1,200);
% 
% hh = LagrangeVal(xx,N,1);
% pphi = hh'*phi*hh;
% 
% subplot(1,2,2)
% pcolor(xx,xx,pphi)
% shading interp
% xlabel('x')
% ylabel('t')
% colorbar
% set(gca,'clim',[-.1 .6])
% axis('square')
% 
% figure
% subplot(1,2,1)
% plot([-1 0 linspace(0,1,100)],[0 0 1/4-1/4*cos(2*pi*linspace(0,1,100))],'--r')
% hold on
% plot(xx,pphi(end,:),xiGLL,phi(end,:),'ob','markerface','b')
% 
% subplot(1,2,2)
% for i=1:200
%     plot(xx,pphi(i,:))%,xiGLL,phi(i,:),'ob','markerface','b')
%     pause(0.1)
% end
