clear all; close all; clc;
N = 5; % nr of cells

% Hyman/Shaskov parameter
c = 0.;

% xi = [ -3/4+GLLnodes(N)/4 -1/4+GLLnodes(N)/4 1/4+GLLnodes(N)/4 3/4+GLLnodes(N)/4 ];
% xi = [ -5/6+GLLnodes(N)/6 -3/6+GLLnodes(N)/6 -1/6+GLLnodes(N)/6 1/6+GLLnodes(N)/6 3/6+GLLnodes(N)/6 5/6+GLLnodes(N)/6];
xi = [ -4/5+GLLnodes(N)/5 -2/5+GLLnodes(N)/5 GLLnodes(N)/5 2/5+GLLnodes(N)/5 4/5+GLLnodes(N)/5];

eta = xi;

% xi_N  = linspace(-1,1,N+1);
% xi_N = GLLnodes(N);
xi_N = xi;
eta_N = xi_N;

x = zeros(5*(N+1)); y = zeros(5*(N+1));
for i=1:5*(N+1)
    for j=1:5*(N+1)
        x(i,j) = xi_N(i)  + c*sin(pi*xi_N(i))*sin(pi*eta_N(j));
        y(i,j) = eta_N(j) + c*sin(pi*xi_N(i))*sin(pi*eta_N(j));
%         x(i,j) = 1/2+1/2*(xi_N(i) + c*cos(pi*xi_N(i)*2).*sin(pi*eta_N(j)*2));
%         y(i,j) = 1/2+1/2*(eta_N(j)+ c*sin(pi*xi_N(i)*2).*cos(pi*eta_N(j)*2));
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
x_xi = zeros(5*(N+1),500); y_xi = x_xi; x_eta = x_xi; y_eta = x_xi;
for i=1:5*(N+1)
        % xi const
        x_xi(i,:)  = xi_N(i)  + c*sin(pi*xi_N(i))*sin(pi*l);
        y_xi(i,:)  = l + c*sin(pi*xi_N(i))*sin(pi*l);
        % eta const
        x_eta(i,:) = l + c*sin(pi*l)*sin(pi*eta_N(i));
        y_eta(i,:) = eta_N(i) + c*sin(pi*l)*sin(pi*eta_N(i));

% % xi const
% x_xi(i,:)  = 1/2+1/2*(xi_N(i) + c*cos(pi*xi_N(i)*2).*sin(pi*l*2));
% y_xi(i,:)  = 1/2+1/2*(l + c*sin(pi*xi_N(i)*2).*cos(pi*l*2));
% % eta const
% x_eta(i,:) = 1/2+1/2*(l + c*cos(pi*l*2).*sin(pi*eta_N(i)*2));
% y_eta(i,:) = 1/2+1/2*(eta_N(i)+ c*sin(pi*l*2).*cos(pi*eta_N(i)*2));

end

% subplot(1,2,2)
hold on
% eta const
for i=1:5*(N+1)
    if rem(i-1,N+1)
        kleur = 'k';
        dikte = 0.5;
    plot(x_xi(i,:),y_xi(i,:),kleur,'linewidth',dikte)
    plot(x_eta(i,:),y_eta(i,:),kleur,'linewidth',dikte)
    end
end
for i=1:5*(N+1)
    if ~rem(i-1,N+1) || i==5*(N+1)
        kleur = 'r';
        dikte = 1;
    plot(x_xi(i,:),y_xi(i,:),kleur,'linewidth',dikte)
    plot(x_eta(i,:),y_eta(i,:),kleur,'linewidth',dikte)
    end
end
axis('square')
% title(['N = ' num2str(N)])
