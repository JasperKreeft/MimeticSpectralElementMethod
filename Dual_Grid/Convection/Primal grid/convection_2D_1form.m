% LIJKT NIET TE WERKEN !!!

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%

% test 1: cosine hill
% test 2: full cosine hill
% test 3: step function
% test 4: cosine function

test = 4;

N = 2;

a = .9;
b = 1;  % MOET 1 BLIJVEN !!!

nr_of_nodes = (N+1)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D  = div(N);
NG = normalgrad(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xiGLL  = GLLnodes(N);
etaGLL = xiGLL;

[h,dhdxi] = LagrangeVal(xiGLL,N,1);
e         = EdgeVal(dhdxi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IT = [ -a*kron(speye(N+1),e') zeros((N+1)^2,N*(N+1))
       zeros((N+1)^2,N*(N+1)) b*kron(speye(N+1),e') ];

C1 = spalloc(nr_of_nodes,2*nr_of_nodes,2*nr_of_nodes);
for j=1:N+1
    C1((j-1)*(N+1)+(1:N+1),(j-1)+(1:N+1:nr_of_nodes)) = speye(N+1);
end
C1(1:nr_of_nodes,nr_of_nodes+(1:nr_of_nodes)) = speye(nr_of_nodes);

C2 = spalloc(N*(N+1),N*(N+1),N*(N+1));
for i=1:N+1
    for j=1:N
        ind1 = j+(i-1)*N;
        ind2 = i+(j-1)*(N+1);
        C2(ind1,ind2) = 1;
    end
end
C2 = kron(speye(2),C2);

X01 = C1*IT*C2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IT2 = zeros(N*(N+1),N^2);
for i=1:N
    IT2((i-1)*(N+1)+(1:N+1),:) = kron(e',[zeros(1,i-1) 1 zeros(1,N-i)]);
end

X12 = [ a*kron(speye(N),e') ; b*IT2 ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_bc = zeros(2*N,1);
for i=1:N

    if test==1
        % cosine hill
        if xiGLL(i)>0
            q_bc(N+i,1) = 0;
        elseif xiGLL(i+1)>0
            q_bc(N+i,1) = -1/4*xiGLL(i)+1/(8*pi)*sin(2*pi*xiGLL(i));
        else
            q_bc(N+i,1) = 1/4*(xiGLL(i+1)-xiGLL(i))-1/(8*pi)*(sin(2*pi*xiGLL(i+1))-sin(2*pi*xiGLL(i)));
        end
    elseif test==2
        % full cosine hill
        q_bc(N+i,1) = 1/4*(xiGLL(i+1)-xiGLL(i))+1/(4*pi)*(sin(pi*xiGLL(i+1))-sin(pi*xiGLL(i)));
        
    elseif test==3
        xstep = -0.4;
        % step-function
        if xiGLL(i+1)<=xstep
            q_bc(N+i,1) = 0;
        elseif xiGLL(i)<=xstep && xiGLL(i+1)>xstep
            q_bc(N+i,1) = xiGLL(i+1)-xstep;
        elseif xiGLL(i)>xstep
            q_bc(N+i,1) = xiGLL(i+1)-xiGLL(i);
        end
    elseif test==4
        % cosine function
        q_bc(N+i,1) = 1/4*(xiGLL(i+1)-xiGLL(i))+1/(8*pi)*(sin(2*pi*xiGLL(i+1))-sin(2*pi*xiGLL(i)));
        if a~=0
        q_bc(i,1) = 1/4*(etaGLL(i+1)-etaGLL(i))-1/(8*pi*a)*(sin(-2*pi*a*(etaGLL(i+1)+1))-sin(-2*pi*a*(etaGLL(i)+1)));
        end
    end
end

% !!!!!
q_bc(N+1:2*N) = -q_bc(N+1:2*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A = X12*D+NG*X01;
% A = X12*D;
% A = NG*X01;

ind_bc = 1:(N+1):N*(N+1);
ind_ic = N*(N+1)+(1:(N+1):N*(N+1));
A([ind_bc ind_ic],:) = [];
A_bc = A(:,[ind_bc ind_ic]);
A(:,[ind_bc ind_ic]) = [];


f = -A_bc*q_bc;

q = A\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qxi = zeros(N+1,N);
Qxi(1,:) = q_bc(1:N);
Qxi(2:N+1,:) = reshape(q(1:N^2),N,N);

Qeta = zeros(N,N+1);
Qeta(:,1) = q_bc(N+(1:N));
Qeta(:,2:N+1) = reshape(q(N^2+(1:N^2)),N,N)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% single_element_volume_postprocessen

% dxideta = diff(xiGLL)'*diff(etaGLL);
% 
% P = reshape(p,N,N);
% 
% % P = P./dxideta;
% 
% for i=1:N
%     for j=1:N
% surf([xiGLL(i:i+1) ;  xiGLL(i:i+1)],[etaGLL(j) etaGLL(j) ; etaGLL(j+1) etaGLL(j+1)],P(i,j)/dxideta(i,j)*ones(2))
% hold on
%     end
% end
% xlabel('\xi')
% ylabel('\eta')
% view([0 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx = linspace(-1,1,200);
XX = xx'*ones(1,200);
YY = XX';

XiGLL  = xiGLL'*ones(1,N+1);
EtaGLL = ones(N+1,1)*etaGLL;

[hh,dhhdxx] = LagrangeVal(xx,N,1);
ee = EdgeVal(dhhdxx);

qq = -ee'*Qeta*hh+hh'*Qxi*ee;

figure
surf(XX,YY,qq)
% pcolor(XX,YY,pp)
shading interp
xlabel('x')
ylabel('t')
colorbar
% set(gca,'clim',[-.1 .6])
% set(gca,'clim',[0 1])
axis('square')

figure
qqbc1 = ee'*q_bc(1:N);
plot(xx,qqbc1)
axis equal
hold on
qqbc2 = ee'*q_bc(N+1:2*N);
plot(xx,qqbc2,'g')
% plot(linspace(-1,1,100),1/4+1/4*cos(2*pi*linspace(-1,1,100)),'r')