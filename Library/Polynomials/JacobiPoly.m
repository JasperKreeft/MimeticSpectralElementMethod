function [P,b1dP,b1] = JacobiPoly(N,alpha,beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [P,b1dP,b1] = JacobiPoly(N,alpha,beta)
% Jacobi polynomials for order N on a domain [-1,1]
% 
% Written by Jasper Kreeft - 2009
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P = zeros(N+1);

P(1,N+1) = 1;                                   % P_0
P(2,[N N+1]) = [(alpha+beta+2)/2 (alpha-beta)/2]; % P_1
for n=1:N-1
    a1 = 2*(n+1)*(n+alpha+beta+1)*(2*n+alpha+beta);
    a2 = (2*n+alpha+beta+1)*(alpha^2-beta^2);
    a3 = (2*n+alpha+beta)*(2*n+alpha+beta+1)*(2*n+alpha+beta+2);
    a4 = 2*(n+alpha)*(n+beta)*(2*n+alpha+beta+2);
    P(n+2,:) =  a2/a1*P(n+1,:)+...
                a3/a1*[P(n+1,2:N+1) 0]+...
               -a4/a1*P(n,:);
end

b1dP = zeros(N+1,N+2);
for n=1:N
    b2_1 = n*(alpha-beta);      % x^0
    b2_2 = -n*(2*n+alpha+beta); % x^1
    b3   = 2*(n+alpha)*(n+beta);% x^0
    b1dP(n+1,:) = b2_1*[0 P(n+1,:)]+...
                b2_2*[P(n+1,:) 0]+...
                b3*[0 P(n,:)];
end

n = 0:N;
b1(n+1,1) = -(2*n+alpha+beta);   % x^2
b1(n+1,2) = 2*n+alpha+beta;      % x^0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To plot the Jacoby polynomials and their first derivative, uncommand the
% lines below

% kk  = 1000;
% xi  = linspace(-1,1,kk);
% psi = zeros(N+1,kk);
% for i=0:N
%     psi(i+1,:) = polyval(P(i+1,:),xi);
% end
% figure
% plot(xi,psi)
% grid
% legend(num2str((0:N)'),4)
% title(['Jacobi polynomials P^',num2str(alpha),'^,^',num2str(beta),'_n,  with n=0...',num2str(N)])
% 
% 
% dpsi = zeros(N+1,kk);
% for i=0:N
%     dpsi(i+1,:) = polyval(b1dP(i+1,:),xi)./(b1(i+1,1).*xi.^2+b1(i+1,2)+1e-15);
% end
% figure
% plot(xi,dpsi)
% grid
% legend(num2str((0:N)'),4)
% title(['Derivative of Jacobi polynomials, dP^',num2str(alpha),'^,^',num2str(beta),'_n/dx,  with n=0...',num2str(N)])