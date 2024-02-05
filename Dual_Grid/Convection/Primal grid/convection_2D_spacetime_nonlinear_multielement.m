% WAARSCHIJNLIJK NOG NIET AF

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 3;

nr_of_nodes = (N+1)^2;

numColumns = 3;
numRows = 1; r=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localnr = reshape(1:nr_of_nodes,N+1,N+1);

globalnr = zeros(numColumns*N+1,numRows*N+1);
% for j=1:numRows
    for c=1:numColumns
        if r==1 && c==1
            globalnr(1:N+1,1:N+1) = localnr;
        elseif r==1 && c>1
            globalnr((c-1)*N+1+(1:N),1:N+1) = (c-1)*N*(N+1)+N+1+reshape(1:N*(N+1),N,N+1);
        end
    end
% end

% flipud(globalnr')


xiGLL = GLLnodes(N);

X = zeros(numColumns*N+1,1);
Xe = linspace(0,1,numColumns+1);
for i=1:numColumns
    X((i-1)*N+(1:N+1)) = (Xe(i+1)+Xe(i))/2+(Xe(i+1)-Xe(i))/2*xiGLL;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Dp,Gd,Cd,NGp,Cp,Dd,Gp] = topology(N);
clear Dp Gd Cd NGp Cp Dd

G = full(kron(speye(numColumns),Gp));
for c=numColumns:-1:2
    ind1 = (c-2)*nr_of_nodes+(N+1:N+1:nr_of_nodes);
    ind2 = (c-1)*nr_of_nodes+(1:N+1:nr_of_nodes);
    G(:,ind1) = G(:,ind1)+G(:,ind2);
    G(:,ind2) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


phi_ic = zeros(numColumns*N+1,1);
for i=1:numColumns*N+1
    phi_ic(i,1) = 1/4-1/4*cos(2*pi*X(i));
end
% plot(X,phi_ic,'-o')


phi_new = zeros((numColumns*N+1)*N+1,1);
for j=1:N+1
    phi_new(globalnr(1:(numColumns*N+1),j),1) = phi_ic;
end

phi0 = ones((numColumns*N+1)*N+1,1);

break

% HIERRRRRRRRRRR

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,dhdxi] = LagrangeVal(xiGLL,N,1);
e         = EdgeVal(dhdxi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IT = [ kron(speye(N+1),e') zeros((N+1)^2,N*(N+1))
       zeros((N+1)^2,N*(N+1)) kron(speye(N+1),e') ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C1 = spalloc(nr_of_nodes,2*nr_of_nodes,2*nr_of_nodes);

C1(1:nr_of_nodes,1:nr_of_nodes) = speye(nr_of_nodes);
for j=1:N+1
    C1((j-1)*(N+1)+(1:N+1),nr_of_nodes+(j-1)+(1:N+1:nr_of_nodes)) = speye(N+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IT = C1*IT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while max(abs(phi_new-phi0))>1e-4

    phi0 = phi_new;

    % Picard linearization
    B = ([phi0*ones(1,N*(N+1)) ones((N+1)^2,N*(N+1))]).*IT*Gp;

    B([1:N+1 (1:N)*(N+1)+1],:) = [];
    f=-B(:,2:N+1)*phi_ic(2:N+1);
    B(:,[1:N+1 (1:N)*(N+1)+1]) = [];

    Phi = B\f;

    phi_new = zeros((N+1)^2,1);
    phi_new(1:N+1,1) = phi_ic;
    phi_new(N+2:N+1:(N+1)^2,1) = zeros(N,1);
    ind = N+2:(N+1)^2; ind(1:N+1:N*(N+1)) = [];
    phi_new(ind) = Phi;

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
% figure
% pause(1)
% for i=1:nn
% %     figure(2)
% %     subplot(1,2,2)
%     plot(xx,pphi(i,:))
%     axis([-1 1 -0.1 0.6])
%     grid on
%     F = getframe;
%     aviobj = addframe(aviobj,F);
% %     pause(0.05)
% end
% plot(xx,pphi(end,:),xiGLL,phi(end,:),'ob','markerface','b')
% axis([-1 1 -0.1 0.6])
% grid on
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
