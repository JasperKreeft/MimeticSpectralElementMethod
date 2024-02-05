clear all
clf%ose all
clc

a = 1;
dt = .02;

N = 10;
T = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL]=GLLnodes(N);

[h,e] = MimeticpolyVal(xiGLL,N,1);

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

G = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N); G(N,1) = 1;

IGLL = e(:,1:N);
IGLL(:,1) = IGLL(:,1)+e(:,N+1);

wp = wGLL(1:N);
wp(1) = wGLL(1)+wGLL(N+1);  % Wel of niet ???
W = spdiags(wp',0,N,N);

A = a*W*IGLL'*G;

% break

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % exact

    u_exact = zeros(1,N);
    for i=1:N
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
    plot(xiGLL,[u_exact u_exact(1)],'-xr')
    hold on
    %%%%%%%%%%%%%%%%%%%%%

    % BE
%     phinew = inv(W+dt*A)*(W*phi);

    % FE
%     phinew = inv(W)*((W-dt*A)*phi);

    % CN
    phinew = inv(W+dt/2*A)*(W-dt/2*A)*phi;

    plot(xiGLL,[phinew ; phinew(1)],'-o')
    ylim([-.2 1])
    pause(0.05)
    hold off
    phi = phinew;
    
end