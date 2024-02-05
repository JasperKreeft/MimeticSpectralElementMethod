N = 6; % nr of cells

xi  = linspace(-1,1,N+1);
% xi = GLLnodes(N);
eta = xi;

% rotation parameters
a = 1/2*sqrt(3);
if abs(a)>1; disp('warning: |a|>1'); end 
b = sqrt(1-a^2);

% Hyman/Shaskov parameter
c = 0.2;

x = zeros(N+1); y = zeros(N+1);
for i=1:N+1
    for j=1:N+1
        x(i,j) = xi(i)  + c*sin(pi*xi(i))*sin(pi*eta(j));
        y(i,j) = eta(j) + c*sin(pi*xi(i)).*sin(pi*eta(j));
%         x(i,j) = xi(i)^3;
%         y(i,j) = eta(j)^3;
%         x(i,j) = xi(i)*sqrt(0.9)  + eta(j)*sqrt(0.1);
%         y(i,j) = eta(j)*sqrt(0.9) + xi(i)*sqrt(0.1);
%         x(i,j) = a*xi(i)  + b*eta(j);
%         y(i,j) = a*eta(j) - b*xi(i);
    end
end


figure
hold on
for i=1:N
    for j=1:N+1
        plot([x(i,j) x(i+1,j)],[y(i,j) y(i+1,j)],'.-')
    end
end
for i=1:N+1
    for j=1:N
        plot([x(i,j) x(i,j+1)],[y(i,j) y(i,j+1)],'.-')
    end
end