clear all
clf%ose all
clc

% if ispc; figure('windowstyle','docked'); else figure; end

v = 1;
dt = .02;

N = 20;
T = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

nn = 1000;
xx = linspace(-1,1,nn);
[hh,ee] = MimeticpolyVal(xx,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = zeros(N,1);
for i=1:N
    if xi(i+1)<=0
        phi(i,1) = 1/4*(xi(i+1)-xi(i))-1/(8*pi)*(sin(2*pi*xi(i+1))-sin(2*pi*xi(i)));
    elseif xi(i+1)>0 && xi(i)<0
        phi(i,1) = 1/4*(0-xi(i))-1/(8*pi)*(0-sin(2*pi*xi(i)));
    else
        phi(i,1) = 0;
    end
end

phidx = phi./diff(xi)';
plot([xi(2:N+1) ; xi(1:N)],[phidx phidx]','g')
hold on
grid on

pphi = phi'*ee;
plot(xx,pphi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N); D(N,1) = 1;
% E = e(:,1:N);
% E(:,1) = E(:,1)+e(:,N+1);
% wp = w(1:N)';
% W = spdiags(wp,0,N,N);
% A = E*W*E'*D*E';
% K2 = E*W*E';
% K1 = M2*D*E';


D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);
E = e;
W = spdiags(w',0,N+1,N+1);
M1 = E*W*E';
K1 = M1*D*E';
K2 = M1;


t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % nodally exact, interpolated solution

    u_interp = zeros(1,N);
    for i=1:N
        s = xi(i)-v*j*dt;
        while s<-1 || s>1
            if s>1
                s = s-2;
            elseif s<0
                s = s+2;
            end
        end

        if s<=0
            u_interp(i) = 1/4-1/4*cos(2*pi*s);
        else
            u_interp(i) = 0;
        end
    end
    
    uu_interp = [u_interp u_interp(1)]*hh;
    
    plot(xi,[u_interp u_interp(1)],'xr')
    hold on
    plot(xx,uu_interp,'r')
    ylim([-1 1])
    %%%%%%%%%%%%%%%%%%%%%

    phinew = timemarching(K2,K1,dt,phi,'ESDIRK3');
    
    pphinew = phinew'*ee;

    phidx = phinew./diff(xi)';
    plot([xi(2:N+1) ; xi(1:N)],[phidx phidx]','g')
    plot(xx,pphinew)
%     ylim([-.2 1])
    pause(0.05)
    hold off
    phi = phinew;
    
end