clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%

numRows    = 1;
numColumns = 2;

N = 12;
N2 = N*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D,G] = topology(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xiG = Gnodes(N);
xiEG = [-1 xiG 1];

[hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);
ew_w           = EdgeVal(dhwdxiw);

IT = kron(speye(N),ew_w');

nr_of_nodes = N2+4*N;
nr_of_lines = 2*N*(N+1);

IT = [ IT zeros(N2+2*N,N*(N+1))
       zeros(N2+2*N,N*(N+1)) IT ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C1 = spalloc(nr_of_nodes,2*(N2+2*N),2*(N2+2*N));

for j=1:N
    ind1 = (j-1)*N+(1:N);
    ind2 = (j-1)*(N+2)+1+(1:N);
    C1(ind1,ind2) = speye(N);
end

ind1 = N2+(1:2:2*N);
ind2 = 1:N+2:N*(N+2);
C1(ind1,ind2) = speye(N);

ind1 = N2+(2:2:2*N);
ind2 = (1:N)*(N+2);
C1(ind1,ind2) = speye(N);

for i=1:N
    ind1 = (1:N:N2)+(i-1);
    ind2 = N*(N+2)+1+(1:N)+(i-1)*(N+2);
    C1(ind1,ind2) = speye(N);
end

ind1 = N2+2*N+(1:2:2*N);
ind2 = N*(N+2)+(1:N+2:N*(N+2));
C1(ind1,ind2) = speye(N);

ind1 = N2+2*N+1+(1:2:2*N);
ind2 = N*(N+2)+(1:N)*(N+2);
C1(ind1,ind2) = speye(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A = C1*IT*G;

IT = C1*IT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_icEG = zeros(N+2,1);
for i=1:N+2
    phi_icEG(i,1) = 1/4+0.00001+1/4*cos(pi*xiEG(i));
end
phi_ic = phi_icEG(2:N+1);

phi_new = zeros(nr_of_nodes,1);
phi_new(1:N2) = kron(ones(N,1),phi_ic);
phi_new(N2+(1:2:2*N),1) = phi_icEG(1); % left
phi_new(N2+(2:2:2*N),1) = phi_icEG(N+2); % right
phi_new(N2+2*N+(1:2:2*N),1) = phi_ic; % lower
phi_new(N2+2*N+(2:2:2*N),1) = phi_ic; % upper

phi0 = ones(nr_of_nodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error = 1;
while error>1e-1

    phi0 = phi_new;

    % Picard linearization
    B = ([phi0*ones(1,N*(N+1)) ones(nr_of_nodes,N*(N+1))]).*IT;
%     break
    B = B*G;

    B([N2+(1:2:2*N) N2+2*N+(1:2:2*N)],:) = [];
    f=-B(:,N2+2*N+(1:2:2*N))*phi_ic;
    B(:,[N2+(1:2:2*N) N2+2*N+(1:2:2*N)]) = [];

    Phi = B\f;

    phi_new = zeros(nr_of_nodes,1);
    phi_new(1:N2) = Phi(1:N2);
    phi_new(N2+(1:2:2*N),1) = zeros(N,1); % left
    phi_new(N2+(2:2:2*N),1) = Phi(N2+(1:N)); % right
    phi_new(N2+2*N+(1:2:2*N),1) = phi_ic; % lower
    phi_new(N2+2*N+(2:2:2*N),1) = Phi(N2+N+(1:N)); % upper

    error = max(abs(phi_new-phi0))

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PHI = zeros(N+2);
PHI(2:N+1,2:N+1) = reshape(Phi(1:N2),N,N);
PHI(2:N+1,1)     = phi_ic;
PHI(1,2:N+1)     = zeros(1,N);
PHI(N+2,2:N+1)   = Phi(N2+(1:N));
PHI(2:N+1,N+2)   = Phi(N2+N+(1:N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nn = 200;
xx = linspace(-1,1,nn); yy = xx;
X = xx'*ones(1,nn);
Y = ones(nn,1)*yy;

hG   = LagrangeVal(xx,N,2);
hEG  = LagrangeVal(xx,N,3);
hGEG = (hG+hEG(2:N+1,:))/2;

pphi = hEG(2:N+1,:)'*PHI(2:N+1,2:N+1)*hEG(2:N+1,:)+... % inner part
       hEG(1,:)'*PHI(1,2:N+1)*hGEG+...                 % left side
       hEG(N+2,:)'*PHI(N+2,2:N+1)*hGEG+...             % right side
       hGEG'*PHI(2:N+1,1)*hEG(1,:)+...                 % lower side
       hGEG'*PHI(2:N+1,N+2)*hEG(N+2,:);                % upper side

% figure
subplot(1,2,1)
surf(X,Y,pphi,'FaceColor','interp',...
	'EdgeColor','none',...
	'FaceLighting','phong')
shading interp
xlabel('x')
ylabel('y')
% colorbar
set(gca,'clim',[-.1 0.6])
% set(gca,'clim',[-.1 1.1])
axis tight
view(-22,50)
camlight left

% figure
subplot(1,2,2)
pcolor(X,Y,pphi)
shading interp
xlabel('x')
ylabel('y')
colorbar
set(gca,'clim',[-.1 .6])
% set(gca,'clim',[-.1 1.1])
axis('square')


% clear aviobj
% aviobj = avifile('example.avi');
% aviobj.KeyFramePerSec = 20;
% aviobj.quality = 20;
figure
% pause(1)
for i=1:nn
%     figure(2)
%     subplot(1,2,2)
    plot(xx,pphi(:,i))
    axis([-1 1 -0.1 0.6])
    grid on
%     F = getframe;
%     aviobj = addframe(aviobj,F);
    pause(0.01)
end
% plot(xx,pphi(end,:),xiGLL,phi(end,:),'ob','markerface','b')
axis([-1 1 -0.1 0.6])
grid on
% F = getframe;
% aviobj = addframe(aviobj,F);
% aviobj = close(aviobj);