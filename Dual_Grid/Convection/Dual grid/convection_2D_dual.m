clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%

% test 1: cosine hill
% test 2: full cosine hill
% test 3: step function
% test 4: cosine function

test = 2;

N = 20;
N2 = N*N;

a = 1;


%%%%%%%%%%%%%%%%%%%%%%

[D,G] = topology(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xiG = Gnodes(N);     etaG  = xiG;
xiEG = [-1 xiG 1];   etaEG = xiEG;

[hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);
ew_w           = EdgeVal(dhwdxiw);

IT = kron(speye(N),ew_w');

nr_of_nodes = N2+4*N;
nr_of_lines = 2*N*(N+1);

IT = [ a*IT zeros(N2+2*N,N*(N+1))
       zeros(N2+2*N,N*(N+1)) IT ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = C1*IT*G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_ic = zeros(N,1);
phi_bc = zeros(N,1);

if test==1
    % cosine hill
    for i=1:N
        if xiG(i)<=0
            phi_ic(i,1) = 1/4-1/4*cos(2*pi*xiG(i));
        else
            phi_ic(i,1) = 0;
        end
    end
elseif test==2
    % full cosine hill
    phi_ic = 1/4+1/4*cos(pi*xiG');
elseif test==3
    % step-function
    xistep = -0.4;
    phi_ic = (xiG>xistep)';
elseif test==4
    % cosine function
    phi_ic = 1/4+1/4*cos(2*pi*xiG');
    phi_bc = 1/4+1/4*cos(-2*pi*a*(etaG'+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B = A;
ind_bc = N2+(1:2:2*N);
ind_ic = N2+2*N+(1:2:2*N);
B([ind_bc ind_ic],:) = [];
f=-B(:,[ind_bc ind_ic])*[phi_bc ; phi_ic];
B(:,[ind_bc ind_ic]) = [];

phi = B\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PHI = zeros(N+2);
PHI(2:N+1,2:N+1) = reshape(phi(1:N2),N,N);
PHI(2:N+1,1)     = phi_ic;
PHI(1,2:N+1)     = phi_bc';
PHI(N+2,2:N+1)   = phi(N2+(1:N));
PHI(2:N+1,N+2)   = phi(N2+N+(1:N));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
