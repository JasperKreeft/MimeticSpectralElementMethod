% JACOBIAN'S AND OTHER TRANSFORMATION MATRICES MUST BE INCLUDED

clear all
close all
clc

if ispc; figure('windowstyle','docked'); else figure; end

N = 5;
H = 6;
v = 1;
T = 10;
L = 1;
dt = .01;

InitInfo.shape = 'cosine3'; %'step'; %
InitInfo.X0    = [0.25 0.75]; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

nn = 1000;
xixi = linspace(-1,1,nn+1);
[hh,ee] = MimeticpolyVal(xixi,N,1);
xx = linspace(0,1,H*nn+1);

dxdxi = L/(2*H);

x = [];
for h=1:H
    xe = (2*h-1)*dxdxi+xi*dxdxi;
    x  = [ x(1:end-1) xe ]; 
end
plot(x,zeros(1,N*H+1),'sk')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(x,InitInfo);

phidx = phi./diff(x)';
plot([x(2:N*H+1) ; x(1:N*H)],[phidx phidx]','g')
hold on
grid on


for h=1:H
    ind = (h-1)*N+(1:N);
    pphi = phi(ind)'*ee/dxdxi;
    hold on
    ind = (h-1)*nn+(1:nn+1);
    plot(xx(ind),pphi,'r')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N*H,1) ones(N*H,1)],[0 1],N*H,N*H); D(N*H,1) = 1;

E1 = zeros(N*H,N*H+1);
for h=1:H
    E1((h-1)*N+(1:N),(h-1)*N+(1:N+1)) = e;
end
% Periodic BC
E1(:,1) = E1(:,1)+E1(:,N*H+1);
E1(:,N*H+1) = [];

% Central
for h=1:H
    E1(:,(h-1)*N+1) = E1(:,(h-1)*N+1)/2;
end
% % Upwind
% for h=1:H-1
%     E1(h*N+(1:N),h*N+1) = zeros(N,1);
% end
% % Downwind
% for h=1:H-1
%     E1((h-1)*N+(1:N),h*N+1) = zeros(N,1);
% end

I1 = zeros(N);
for i=1:N
    for j=1:N
        I1(i,j) = sum(w.*e(i,:).*e(j,:));
    end
end
I1 = kron(eye(H),I1); % Only possible for uniform grid !!!

A = I1*D*E1'/dxdxi;

% Ebc = zeros(N*H,N*H+1);
% Ebc(1:N,1) = e(:,1);
% Ebc(N*(H-1)+(1:N),N*H+1) = e(:,N+1);
% B = Ebc*Ebc'/dxdxi;

M1 = v*0.5*(A-A');
% M1 = v*A;
M2 = I1;


t=0; j=0;
ref.shape = 'cosine2'; ref.X0 = InitInfo.X0;
while t<T

    t = t+dt;
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % nodally exact, interpolated solution

    u_interp = zeros(1,N*H);
    for i=1:N*H
        s = x(i)-v*j*dt;
        while s<0 || s>L
            if s>L
                s = s-L;
            elseif s<0
                s = s+L;
            end
        end

        u_interp(i) = initialfunction(s,ref);
    end

    u_interp = [u_interp u_interp(1)];
    for h=1:H
        ind = (h-1)*N+(1:N+1);
        uu_interp = u_interp(ind)*hh;
    
        ind = (h-1)*nn+(1:nn+1);
        plot(xx(ind),uu_interp,'r','linewidth',4)
        hold on
    end
   
    %%%%%%%%%%%%%%%%%%%%%

    phinew = timemarching(M2,M1,dt,phi,'ESDIRK3');
    
%     if ~exist('phiold','var')
%         phinew = timemarching(M2,M1,dt,phi,'BE');
%     else
%         phinew = timemarching(M2,M1,dt,[phi phiold],'BDF2');
%     end

    for h=1:H
        ind1 = (h-1)*N+(1:N);
        ind2 = (h-1)*N+(1:N+1);
        pphinew = phinew(ind1)'*ee/dxdxi;

        phidx = phinew(ind1)./diff(x(ind2))';
        plot([x(ind1+1) ; x(ind1)],[phidx phidx]','g','linewidth',4)
        ind = (h-1)*nn+(1:nn+1);
        plot(xx(ind),pphinew)
        ylim([-.2 1.2])
    end
    pause(0.15)
    hold off
    
    phiold = phi;
    phi    = phinew;

end