function [F] = force(x,y)

global N problem c

n = 8;

[xibar,Gw] = Gnodes(n); GGw = Gw'*Gw;
etabar = xibar;
F = zeros(N);
for i=1:N
    for j=1:N
        xi  = ( (x(i)+x(i+1))/2+(x(i+1)-x(i))/2*xibar )'*ones(1,n);
        eta = ones(n,1)*( (y(j)+y(j+1))/2+(y(j+1)-y(j))/2*etabar );

        xx = xi +c*sin(pi*xi).*sin(pi*eta);
        yy = eta+c*sin(pi*xi).*sin(pi*eta);
        
        k11 = ones(n);
        k12 = zeros(n);
        k21 = zeros(n);
        k22 = ones(n);
        
        dk11 = zeros(n);
        dk12 = zeros(n);
        dk21 = zeros(n);
        dk22 = zeros(n);

        switch problem
            case 'sine'
                global m                                         %#ok<TLEV>
                f = -m^2*pi^2*(k11+k22).*sin(m*pi*xx).*sin(m*pi*yy)+...
                    +m^2*pi^2*(k12+k21).*cos(m*pi*xx).*cos(m*pi*yy)+...
                    +m*pi*(dk11+dk21).*cos(m*pi*xx).*sin(m*pi*yy)+...
                    +m*pi*(dk12+dk22).*sin(m*pi*xx).*cos(m*pi*yy);
            case 'cosine'
                f = -pi^2/4*(cos(pi*xi).*(cos(pi*eta)+1)+(cos(pi*xi)+1).*cos(pi*eta));
        end

        dxdxi  = 1+c*pi*cos(pi*xi).*sin(pi*eta);
        dxdeta = c*pi*sin(pi*xi).*cos(pi*eta);
        dydxi  = c*pi*cos(pi*xi).*sin(pi*eta);
        dydeta = 1+c*pi*sin(pi*xi).*cos(pi*eta);

        J = dxdxi.*dydeta-dxdeta.*dydxi;
        
        Jbar = ((x(i+1)-x(i))/2)*((y(j+1)-y(j))/2)*ones(n);

        F(i,j) = sum(sum(GGw.*f.*J.*Jbar));
        
    end
end

% sum(sum(F))