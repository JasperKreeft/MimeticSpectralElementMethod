clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%

N = 20;
M = N;

%%%%%%%%%%%%%%%%%%%%%%

[Dp,Gd,Cd,NGp,Cp,Dd,Gp] = topology(N);
clear Dp Gd Cd NGp Cp Dd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xiGLL = GLLnodes(N);

[h,dhdxi] = LagrangeVal(xiGLL,N,1);
e         = EdgeVal(dhdxi);

phi_ic = zeros(N+1,1);
for i=1:N+1
%     if xiGLL(i)<=0
%         phi_ic(i,1) = 1/4-1/4*cos(2*pi*xiGLL(i));
%     else
%         phi_ic(i,1) = 0;
%     end

    phi_ic(i,1) = 1/4+1/4*cos(pi*xiGLL(i));
    
%     if xiGLL(i)>-0.4
%         phi_ic(i,1) = 1;
%     end
end

phi_new =  kron(ones(N+1,1),phi_ic);
phi0 = ones((N+1)^2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = kron(speye(N+1),e);
IT = I';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_of_nodes = (N+1)^2;
C1 = spalloc(nr_of_nodes,nr_of_nodes,nr_of_nodes);

for j=1:N+1
    C1((j-1)*(N+1)+(1:N+1),(j-1)+(1:N+1:nr_of_nodes)) = speye(N+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gx = Gp(1:N*(N+1),:);
Gt = Gp(N*(N+1)+1:2*N*(N+1),:);

B1 = IT*Gx;
IGt = (C1*IT)*Gt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error = 1;
while error>1e-4

    phi0 = phi_new;

    % Newton linearization !!!! WERKT NOG NIET GOED !!!!
    A = (-Gx.*(ones(N*(N+1),1)*phi0'))'*I+...
        (I.*(ones(N*(N+1),1)*phi0'))'*Gx+... 
        IGt;

    f1 = ((phi0*ones(1,N*(N+1))).*IT)*(Gx*phi0);

    A([1:N+1 (1:N)*(N+1)+1],:) = [];
    f1([1:N+1 (1:N)*(N+1)+1])  = [];

    f = f1-A(:,2:N+1)*phi_ic(2:N+1);

    A(:,[1:N+1 (1:N)*(N+1)+1]) = [];

    Phi = A\f;

    phi_new = zeros((N+1)^2,1);
    phi_new(1:N+1,1) = phi_ic;
    phi_new(N+2:N+1:(N+1)^2,1) = zeros(N,1);
    ind = N+2:(N+1)^2; ind(1:N+1:N*(N+1)) = [];
    phi_new(ind) = Phi;

    error = max(abs(phi_new-phi0))

    if error > 10
        break
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = reshape(Phi,N,N)';

phi = [phi_ic(2:N+1)' ; phi ];

phi = [zeros(N+1,1) phi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(1,2,1)
% pcolor(xiGLL,xiGLL,phi)
% shading interp
% xlabel('x')
% ylabel('t')
% colorbar
% set(gca,'clim',[-.1 .6])
% % set(gca,'clim',[-.1 1.1])
% axis('square')

nn = 400;
xx = linspace(-1,1,nn);

hh = LagrangeVal(xx,N,1);
pphi = hh'*phi*hh;

% subplot(1,2,2)
pcolor(xx,xx,pphi)
shading interp
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 .6])
% set(gca,'clim',[-.1 1.1])
axis('square')

figure
% subplot(1,2,1)
plot([-1 0 linspace(0,1,100)],[0 0 1/4-1/4*cos(2*pi*linspace(0,1,100))],'--r')
hold on
plot(xx,pphi(end,:),xiGLL,phi(end,:),'ob','markerface','b')

% clear aviobj
% aviobj = avifile('example.avi');
% aviobj.KeyFramePerSec = 20;
% % aviobj.quality = 20;
figure
pause(1)
for i=1:nn
%     figure(2)
%     subplot(1,2,2)
    plot(xx,pphi(i,:))
    axis([-1 1 -0.1 0.6])
    grid on
%     F = getframe;
%     aviobj = addframe(aviobj,F);
    pause(0.02)
end
plot(xx,pphi(end,:),xiGLL,phi(end,:),'ob','markerface','b')
axis([-1 1 -0.1 0.6])
grid on
% F = getframe;
% aviobj = addframe(aviobj,F);
% aviobj = close(aviobj);

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
