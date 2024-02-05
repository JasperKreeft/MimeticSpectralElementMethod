clear all
clf%ose all
clc

a = 1;
dt = .01;

N = 30;
T = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL]=GLLnodes(N);
dx = diff(xiGLL);


[h,dhdx] = LagrangeVal(xiGLL,N,1);
e = EdgeVal(dhdx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = zeros(N,1);
for i=1:N
    if xiGLL(i)<=0
        phi(i,1) = 1/4-1/4*cos(2*pi*xiGLL(i));
    else
        phi(i,1) = 0;
    end
end

plot(xiGLL,[phi ; phi(1)],'-o')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N);
G(N,1) = 1;

IGLL = e(:,1:N);
IGLL(:,1) = IGLL(:,1)+e(:,N+1);


A = IGLL'*G;

% break

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % exact

    for i=1:N
        s = xiGLL(i)-a*j*dt;
        while s<-1 || s>1+dx(i)
            if s>1+dx(i)
                s = s-(2+dx(i));
            elseif s<0
                s = s+(2+dx(i));
            end
        end

        if s<=0
            u_exact(i) = 1/4-1/4*cos(2*pi*s);
        else
            u_exact(i) = 0;
        end
    end

    plot(xiGLL,[u_exact u_exact(1)],'-xr')
    hold on
    %%%%%%%%%%%%%%%%%%%%%

    % BE
%     phinew = inv(speye(N)+dt*A)*phi;

    % FE
%     phinew = (speye(N)-dt*A)*phi;

    % CN
    phinew = inv(speye(N)+dt/2*A)*(speye(N)-dt/2*A)*phi;

    plot(xiGLL,[phinew ; phinew(1)],'-o')
    ylim([-.2 1])
    pause(0.1)
    hold off
    phi = phinew;
    
end