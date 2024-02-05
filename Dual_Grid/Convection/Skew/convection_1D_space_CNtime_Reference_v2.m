clear all
clf%ose all
clc

a = 2;
dt = .01;

N = 20;
T = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL]=GLLnodes(N);
dx = diff(xiGLL);


[h,dhdx] = LagrangeVal(xiGLL,N,1);
e = EdgeVal(dhdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = zeros(N+1,1);
for i=1:N+1
    if xiGLL(i)<=0
        phi(i,1) = 1/4-1/4*cos(2*pi*xiGLL(i));
    else
        phi(i,1) = 0;
    end
end

plot(xiGLL,phi,'-o')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);

IGLL = e;

W = spdiags(wGLL',0,N+1,N+1);

A = a*W*IGLL'*G;

% break

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % exact

    u_exact = zeros(1,N+1);
    for i=1:N+1
        s = xiGLL(i)-a*j*dt;
        while s<-1 || s>1
            if s>1
                s = s-2;
            elseif s<0
                s = s+2;
            end
        end

        if s<=0
            u_exact(i) = 1/4-1/4*cos(2*pi*s);
        else
            u_exact(i) = 0;
        end
    end
    plot(xiGLL,u_exact,'-xr')
    hold on
    %%%%%%%%%%%%%%%%%%%%%

    % BE
    phinew = inv(W+dt*A)*(W*phi);

    % FE
%     phinew = (speye(N)-dt*A)*phi;

    % CN
%     phinew = inv(speye(N)+dt/2*A)*(speye(N)-dt/2*A)*phi;

    plot(xiGLL,phinew,'-o')
    ylim([-.2 1])
    pause(0.05)
    hold off
    phi = phinew;
    
end