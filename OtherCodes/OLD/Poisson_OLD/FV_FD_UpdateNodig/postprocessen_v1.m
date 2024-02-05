%% Exact
phi_exact = zeros((N+2));
for j=1:N+2
    for i=1:N+2
        phi_exact(i,j) = sin(m*pi*xd_ex(i))*sin(m*pi*yd_ex(j));
%         phi_exact(i,j) = 1/4*(cos(pi*xd_ex(i))+1)*(cos(pi*yd_ex(j))+1);
    end
end

PHI_EXACT = reshape(phi_exact,(N+2)*(N+2),1);

x = linspace(-1,1,100);
y = linspace(-1,1,100);

xx = x'*ones(1,100);
yy = ones(100,1)*y;

phi_ex = sin(m*pi*xx).*sin(m*pi*yy);
% phi_ex = 1/4*(cos(pi*xx)+1).*(cos(pi*yy)+1);


% surf(xx,yy,phi_ex')
contour(xx,yy,phi_ex','k')
hold on


%%

phi = zeros(N+2); l=0;
for i=1:N+2
    for j=1:N+2
        if i==1 || i==N+2
            phi(i,j) = phi_exact(i,j);
        elseif j==1 || j==N+2
            phi(i,j) = phi_exact(i,j);
        else            
        l=l+1;
        phi(i,j) = phi_in(l);
        end
    end
end

% contour(xd_ex,yd_ex,phi')
% hold off

PHI = reshape(phi,[],1);

%% reconstruction

n = 11;
P = zeros((N+1)*n);
xbb = zeros(1,(N+1)*n);
for i=1:N+1
    xb = linspace(0,xd_ex(i+1)-xd_ex(i),n)/(xd_ex(i+1)-xd_ex(i)); yb = xb;
    xbb(1,(i-1)*n+(1:n)) = xd_ex(i)+(xd_ex(i+1)-xd_ex(i))*xb;
    for j=1:N+1
        P((i-1)*n+(1:n),(j-1)*n+(1:n)) = phi(i,j)+(phi(i+1,j)-phi(i,j))*xb'*ones(1,n)+(phi(i,j+1)-phi(i,j))*ones(n,1)*yb+...
            ((phi(i+1,j+1)-phi(i,j))-(phi(i+1,j)-phi(i,j))-(phi(i,j+1)-phi(i,j)))*xb'*yb;
    end
end
ybb = xbb;

% figure
% surf(xbb,ybb,P)
contour(xbb,ybb,P); colorbar; set(gca,'clim',[-1 1]);
% hold on
plot(xd_ex,-ones(1,N+2),'xg')
plot(xd_ex, ones(1,N+2),'xg')
plot(-ones(1,N),yd,'xg')
plot( ones(1,N),yd,'xg')
mark = 'x'; if N>=10; mark = '.'; end;
for i=1:N
    plot(xd,yd(i)*ones(1,N),[mark 'k'])
end
axis('square')
title(['N = ' num2str(N)])
hold off
pause


%% error

% errorL1(N+1-Z(1)) = 1/N^2*sum(sum(abs(phi-phi_exact)));
% errorL2(N+1-Z(1)) = sqrt(1/N^2*sum(sum((phi-phi_exact).^2)));

% errorL1(N+1-Z(1)) = sum(sum( abs(P-phi_ex)*(2/100)^2 ));
% errorL2(N+1-Z(1)) = sqrt( sum(sum( abs(P-phi_ex).^2*(2/100)^2 )) );
% 
% P_interp = zeros(100);
% 
% errorL2_interp(N+1-Z(1)) = sqrt( sum(sum( abs(P_interp-phi_ex).^2*(2/100)^2 )) );