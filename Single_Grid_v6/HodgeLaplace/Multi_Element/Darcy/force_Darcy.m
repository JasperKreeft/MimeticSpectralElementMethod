function [F] = force_Darcy(Mesh)

global N

n = 2*N;

[xi,Gw] = Gnodes(n); GGw = Gw'*Gw;
eta = xi;

Nx = length(x)-1;
Ny = length(y)-1;

F = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        xx  = ( (x(i)+x(i+1))/2+(x(i+1)-x(i))/2*xi )'*ones(1,n);
        yy = ones(n,1)*( (y(j)+y(j+1))/2+(y(j+1)-y(j))/2*eta );
        
        k11 = 10*ones(n);
        k12 =  4*ones(n);
        k21 =  4*ones(n);
        k22 =    ones(n);
        
        dk11 = zeros(n);
        dk12 = zeros(n);
        dk21 = zeros(n);
        dk22 = zeros(n);

        f = -pi^2*(k11+k22).*sin(pi*xx).*sin(pi*yy)+...
            +pi^2*(k12+k21).*cos(pi*xx).*cos(pi*yy)+...
            +pi*(dk11+dk21).*cos(pi*xx).*sin(pi*yy)+...
            +pi*(dk12+dk22).*sin(pi*xx).*cos(pi*yy);
        
        alpha = 1e-3;
        if x(i)<0
            f = alpha/8*(yy.^2-1+(xx+1).^2);
        else
            f = 1/8*(1-xx.^2+1-yy.^2)+alpha/4*(1-yy.^2+xx.*(1-xx));
        end
      
        J = ((x(i+1)-x(i))/2)*((y(j+1)-y(j))/2)*ones(n);

        F(i,j) = sum(sum(GGw.*f.*J));
        
    end
end

% sum(sum(F))

% surf(F)
% figure