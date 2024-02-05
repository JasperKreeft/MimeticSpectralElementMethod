function [] = PostProcessen(Ne,xg,xe,xen,Pe,c,u,u_exact)

L = xg(end);


% exact solution
plot(u_exact(1,:),u_exact(2,:),'m','linewidth',2)

hold on

% Approximate solution
kk = 100;

U = zeros(Ne,kk);

for i=1:Ne
    p = Pe(i);
    La = LagrangePoly(xen(i,1:p+1));
    x(i,:) = linspace(xe(i),xe(i+1),kk);
    for k=1:kk
        for j=1:p+1
            U(i,k) = U(i,k) + u(c(i,j))*La(j,1)*prod(x(i,k)-La(j,2:p+1));
        end
    end
end
plot(x',U','linewidth',2)

plot(xg, u,'+k','linewidth',2)

plot(xe,[U(:,1); U(end,end)]','sk','linewidth',2)

% Grid points
ye= zeros(Ne+1,1);
plot(xe,ye,'d','linewidth',2)


xlabel('x'); ylabel('u');
xlim([0 L])



