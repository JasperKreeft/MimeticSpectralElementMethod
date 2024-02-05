clear all
close all
clc

a = 1;
dt = .05;

N = 20;
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
    if xiGLL(i+1)<=0
        phi(i,1) = 1/4-1/4*cos(2*pi*xiGLL(i+1));
    else
        phi(i,1) = 0;
    end
end

plot(xiGLL,[0 ; phi],'-o')
pause(0.1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);

IGLLT = e';

A = IGLLT*G;

A = A(2:N+1,2:N+1);


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
    
    plot(xiGLL,[u_exact u_exact(1)],'xr')
    hold on
    %%%%%%%%%%%%%%%%%%%%%

    % BE
%     phinew = inv(speye(N)+dt*A)*phi;

    % FE
%     phinew = (speye(N)-dt*A)*phi;

    % CN
    phinew = inv(speye(N)+dt/2*A)*(speye(N)-dt/2*A)*phi;

    plot(xiGLL,[0 ; phinew],'-o')
    ylim([0 1])
    pause%(0.1)
    hold off
    phi = phinew;
    legend('Exact','Numerical')
    
end