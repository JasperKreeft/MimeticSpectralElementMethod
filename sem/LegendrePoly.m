function [Le,dLe] = LegendrePoly(N)
% Legendre polynomials for order N on a domain [-1,1]
% 
% Written by Jasper Kreeft - 2009
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Le = zeros(N+1);

Le(1,N+1) = 1;
Le(2,[N N+1]) = [1 0];
a = N:-1:1;
for i=3:N+1
    Le(i,N+1) = -(i-2)/(i-1)*Le(i-2,N+1);
    Le(i,a) = (2*i-3)/(i-1)*Le(i-1,a+1)-(i-2)/(i-1)*Le(i-2,a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To plot the Legendre polynomials, uncommand the lines below

% kk = 100;
% x = linspace(-1,1,kk);
% f = zeros(N+1,kk);
% for i=1:N+1
%     f(i,:) = polyval(Le(i,:),x);
% end
% 
% plot(x,f)
% grid
% legend(num2str((0:N)'),4)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Differentiating Le
dLe = zeros(N+1,N);
for i=1:N+1
    dLe(i,:) = Le(i,1:N).*(N:-1:1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Lagrange polynomials using Legendre polynomials

% xj = zeros(N+1);
% xj(1,1) = 1;
% for n=1:N
%     xj(n+1,1:n+1) = gausslobattolegendre(n);
% end
% 
% Le_j = zeros(N+1);
% for n=1:N+1
%     Le_j(n,1:n) = polyval(Le(n,:),xj(n,1:n));
% end
% 
% kk = 100;
% x = linspace(-1,1,kk);
% hp = zeros(N+1,kk);
% for i=1:N+1
%     for k=1:kk
%         hp(i,k) = ((x(k)-1)*(x(k)+1)*polyval(dLe(N+1,:),x(k))) / ((N+1)*N*Le_j(N+1,i)*(x(k)-xj(N+1,i)));
%     end
% end
% figure
% plot(x,hp)
