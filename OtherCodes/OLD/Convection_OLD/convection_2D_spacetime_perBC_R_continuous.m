clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Right boundary is periodic boundary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 30;

% for N=10:5:50

a = 2;

%%%%%%%%%%%%%%%%%%%%%%

[Dp,Gd,Cd,NGp,Cp,Dd,Gp] = topology(N);
clear Dp Gd Cd NGp Cp Dd
% per bc.
Gp(:,1:(N+1):(N+1)^2) = Gp(:,1:(N+1):(N+1)^2)+Gp(:,(N+1):(N+1):(N+1)^2);
Gp(:,(N+1):(N+1):(N+1)^2) = [];


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
C1(1:(N+1):(N+1)^2,:) = C1(1:(N+1):(N+1)^2,:)+C1((N+1):(N+1):(N+1)^2,:);
C1((N+1):(N+1):(N+1)^2,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
A = C1*IT*Gp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


phi_ic = zeros(1,N);
for i=1:N
    if xiGLL(i)>=-0.75 && xiGLL(i)<=0.25
        phi_ic(1,i) = 1/2-1/2*sin(2*pi*xiGLL(i));
    else
        phi_ic(1,i) = 0;
    end

%         if xiGLL(i)>-0.4
%             phi_ic(1,i) = 1;
%         end
end

T = 1000;
t = 0;
dt = 2;
while t<T
    t = t+dt

    B = A;
    B(1:N,:) = [];
    f=-B(:,1:N)*phi_ic';
    B(:,1:N) = [];

    Phi = B\f;

    phi = reshape(Phi,N,N)';

    phi = [phi_ic ; phi ];

    phi_ic = phi(end,:);

    phi = [phi phi(:,1)];

    % figure(1)
    % subplot(1,2,1)
    % pcolor(xiGLL,xiGLL,phi)
    % shading interp
    % xlabel('x')
    % ylabel('t')
    % colorbar
    % set(gca,'clim',[-.1 .6])
    % axis('square')



    xx = linspace(-1,1,200);

    hh = LagrangeVal(xx,N,1);
    pphi = hh'*phi*hh;
    pphi_ex = [zeros(1,25) 1/2-1/2*sin(2*pi*xx(26:125)) zeros(1,75) ];

%     figure(1)
%     surf(xx,xx+(t-1),pphi)
%     view([0 0 1])
%     hold on

 
    
%     pause
end
% figure(1)
% shading interp
% xlabel('x')
% ylabel('t')
% colorbar
% set(gca,'clim',[-.1 .6])
% axis('square')

% 
   figure(2)
    % subplot(1,2,1)
    plot(xx,pphi_ex,'--r')
    hold on
    plot(xx,pphi(end,:),'g')%,xiGLL,phi(end,:),'ob','markerface','b','markersize',8)
    hold off

% for i=1:200
%     figure(2)
%     subplot(1,2,2)
%     plot(xx,pphi(i,:))%,xiGLL,phi(i,:),'ob','markerface','b')
%     pause(0.1)
% end


% figure
% surf(xx,xx,pphi,'FaceColor','interp',...
% 	'EdgeColor','none',...
% 	'FaceLighting','phong')
% shading interp
% xlabel('x')
% ylabel('t')
% colorbar
% set(gca,'clim',[-.1 0.6])
% % set(gca,'clim',[-.1 1.1])
% axis tight
% view(-50,30)
% camlight left

% error(N) = sqrt(sum((pphi(200,:)-pphi_ex).^2));
% 
% 
% end
% error=error(10:5:50)'