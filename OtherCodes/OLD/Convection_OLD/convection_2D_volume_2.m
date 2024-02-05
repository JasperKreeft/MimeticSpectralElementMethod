clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%

% test 1: cosine hill
% test 2: full cosine hill
% test 3: step function
% test 4: cosine function

test = 4;

N = 21;
N2 = N*N;

a = 1;
b = 1;  % MOET 1 BLIJVEN !!!

%%%%%%%%%%%%%%%%%%%%%%

Dp = topology(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xiGLL  = GLLnodes(N);
etaGLL = xiGLL;

[h,dhdxi] = LagrangeVal(xiGLL,N,1);
e         = EdgeVal(dhdxi);

IT2 = zeros(N*(N+1),N^2);
for i=1:N
    IT2((i-1)*(N+1)+(1:N+1),:) = kron(e',[zeros(1,i-1) 1 zeros(1,N-i)]);
end

T3e = spalloc(2*N*(N+1),4*N,4*N);
for i=1:N
    T3e((i-1)*(N+1)+1,2*i-1) = 1;
    T3e(i*(N+1),2*i) = 1;
    
    T3e((i-1+N)*(N+1)+1,2*(i+N)-1) = 1;
    T3e((i+N)*(N+1),2*(i+N)) = 1;
end

X12 = [ [ a*kron(speye(N),e') ; b*IT2 ] T3e ];

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
        q_bc(i,1) = 1/4*(etaGLL(i+1)-etaGLL(i))-1/(8*pi*a)*(sin(-2*pi*a*(etaGLL(i+1)+1))-sin(-2*pi*a*(etaGLL(i)+1)));
    end
end

ind_bc1 = [ 1:(N+1):N*(N+1) N*(N+1)+(1:(N+1):N*(N+1)) ];

ind_bc2 = [ N2+(1:2:2*N) N2+2*N+(1:2:2*N) ];


keyboard

D_bc = Dp(:,ind_bc1);

Dp(:,ind_bc1)  = [];
X12(ind_bc1,:) = [];
A = Dp*X12;
A(ind_bc2,:) = [];
A(:,ind_bc2) = [];
D_bc(ind_bc2,:) = [];


f = -D_bc*q_bc;

p = A\f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% single_element_volume_postprocessen

dxideta = diff(xiGLL)'*diff(etaGLL);

P = reshape(p(1:N2),N,N);

% P = P./dxideta;

for i=1:N
    for j=1:N
surf([xiGLL(i:i+1) ; xiGLL(i:i+1)],[etaGLL(j) etaGLL(j) ; etaGLL(j+1) etaGLL(j+1)],P(i,j)/dxideta(i,j)*ones(2))
hold on
    end
end
xlabel('\xi')
ylabel('\eta')
view([0 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xx = linspace(-1,1,200);
XX = xx'*ones(1,200);
YY = XX';

XiGLL  = xiGLL'*ones(1,N+1);
EtaGLL = ones(N+1,1)*etaGLL;

[hh,dhhdxx] = LagrangeVal(xx,N,1);
ee = EdgeVal(dhhdxx);
pp = ee'*P*ee;

figure
surf(XX,YY,pp)
% pcolor(XX,YY,pp)
shading interp
xlabel('x')
ylabel('t')
colorbar
set(gca,'clim',[-.1 .6])
% set(gca,'clim',[-.1 1.1])
axis('square')

figure
qq = ee'*q_bc(1:N);
plot(xx,qq)
axis equal
hold on
qq = ee'*q_bc(N+1:2*N);
plot(xx,qq,'g')
% plot(linspace(-1,1,100),1/4+1/4*cos(2*pi*linspace(-1,1,100)),'r')