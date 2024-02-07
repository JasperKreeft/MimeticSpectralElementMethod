function [x,w]=EGLnodes(N)

[x,w] = GLnodes(N);
x = [ -1 x 1 ];
w = [ 0 w 0 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% La = LagrangePoly(x,2);
% kk = 1000;
% xx = linspace(min(x),max(x),kk);
% f = zeros(N,kk);
% for i=1:N+2
%     xj    = x;
%     xj(i) = [];
%     for k=1:kk
%         f(i,k) = La(i,1)*prod(xx(k)-xj);
%     end
% end
% 
% figure
% plot(xx,f)
% hold on
% plot(x,zeros(N+2),'o','MarkerFaceColor','k','MarkerEdgeColor','k')
% legend(num2str((1:N+2)'))
% grid on
% xlim([-1 1])
