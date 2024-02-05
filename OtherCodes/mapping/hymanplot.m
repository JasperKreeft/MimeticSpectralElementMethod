clear all; close all; clc;
N = 4; % nr of cells

% Hyman/Shaskov parameter
c = 0.2;

xi = GLLnodes(N);%linspace(-1,1,N+1);%
eta = xi;

% xi_N  = linspace(-1,1,N+1);
xi_N = GLLnodes(N);
eta_N = xi_N;

x = zeros(N+1); y = zeros(N+1);
for i=1:N+1
    for j=1:N+1
        x(i,j) = xi_N(i)  + c*sin(pi*xi_N(i))*sin(pi*eta_N(j));
        y(i,j) = eta_N(j) + c*sin(pi*xi_N(i))*sin(pi*eta_N(j));
    end
end

figure
% subplot(1,2,1)
% hold on
% % eta const
% for i=1:N
%     for j=1:N+1
%         plot([x(i,j) x(i+1,j)],[y(i,j) y(i+1,j)],'.-')
%     end
% end
% % xi const
% for i=1:N+1
%     for j=1:N
%         plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'.-')
%     end
% end
% axis('square')
% title(['N = ' num2str(N)])


l = linspace(-1,1,500);
x_xi = zeros(N+1,500); y_xi = x_xi; x_eta = x_xi; y_eta = x_xi;
for i=1:N+1
        % xi const
        x_xi(i,:)  = xi_N(i)  + c*sin(pi*xi_N(i))*sin(pi*l);
        y_xi(i,:)  = l + c*sin(pi*xi_N(i))*sin(pi*l);
        % eta const
        x_eta(i,:) = l + c*sin(pi*l)*sin(pi*eta_N(i));
        y_eta(i,:) = eta_N(i) + c*sin(pi*l)*sin(pi*eta_N(i));
end

% subplot(1,2,2)
hold on
% eta const
for i=1:N+1
    plot(x_xi(i,:),y_xi(i,:),'k')
    plot(x_eta(i,:),y_eta(i,:),'k')
end
axis('square')
% title(['N = ' num2str(N)])
