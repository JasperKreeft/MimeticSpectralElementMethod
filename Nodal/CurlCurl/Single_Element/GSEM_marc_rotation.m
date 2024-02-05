clc;
clear all;
close all;

% This program solves the curl-curl eigenvalue problem on a rectangular
% domain using nodal spectral elements. The puropse of this code is to
% compare the results with the mimetic spectral element approach.
%
% 19 August 2010, Marc Gerritsma

prob = 1
angle = 53;

N=10;        % Polynomial degree
L=pi;        % Size of the domain
c = 0.0;     % Curvature of the gridlines


[xgl wgl] = GLLnodes(N) ;  % GL-points and weights

%% Legendre polynomial of degree N evaluated at the Gauss-Lobatto points

for i=1:N+1
    Lprev = 1;
    Lthis = xgl(i);
    for k=1:N-1;
        LN(i) = ( (2*k+1)*xgl(i)*Lthis - k*Lprev )/(k+1);
        Lprev = Lthis;
        Lthis = LN(i);
    end
end

% Setting up the derivative matrix

for i=1:N+1
    for j=1:N+1
        if i==j
            if i==1
                dh(i,j) = -N*(N+1)/4;
            elseif i==N+1
                dh(i,j) = N*(N+1)/4;
            else
                dh(i,j) = 0;
            end
        else
            dh(i,j) = LN(j)/( LN(i)*(xgl(j)-xgl(i)));
        end
    end
end

%  Loop over the basis functions
A = zeros((N+1)*(N+1),2*(N+1)*(N+1));
for i=1:N+1
    for j=1:N+1
        % Loop over the integration points
        for p=1:N+1
            for q=1:N+1
                if prob == 1
                    dxdxi  = (L/2)*(1 + L*c*cos(L*xgl(p))*sin(L*xgl(q)));
                    dxdeta = (L/2)*L*c*sin(L*xgl(p))*cos(L*xgl(q));
                    dydxi  = (L/2)*L*c*cos(L*xgl(p))*sin(L*xgl(q));
                    dydeta = (L/2)*(1 + L*c*sin(L*xgl(p))*cos(L*xgl(q)));
                elseif prob == 2
                    dxdxi  = (L/2)*cos(angle);
                    dxdeta = (L/2)*sin(angle);
                    dydxi  = -(L/2)*sin(angle);
                    dydeta = (L/2)*cos(angle);
                end
                J = dxdxi*dydeta - dxdeta*dydxi;
                eqnr = p + (q-1)*(N+1);
                if q==j
                    var  = i + (q-1)*(N+1);
                    A(eqnr,var) = A(eqnr,var) + dydeta*dh(i,p)/J;  % dv/dx :  dv/dxi*dxi/dx = dx/dxi * dy/deta/J
                    var = (N+1)*(N+1) + i + (q-1)*(N+1);
                    A(eqnr,var) = A(eqnr,var) + dxdeta*dh(i,p)/J;
                end
                if p==i
                    var = (N+1)*(N+1) + p + (j-1)*(N+1);
                    A(eqnr,var) = A(eqnr,var) - dxdxi*dh(j,q)/J;
                    var = p + (j-1)*(N+1);
                    A(eqnr,var) = A(eqnr,var) - dydxi*dh(j,q)/J;
                end
            end
        end
    end
end

for i=1:N+1
    for j=1:N+1
        if prob == 1
            vecW(i+(j-1)*(N+1)) = (L^2/4)*(1 + L*c*cos(L*xgl(i))*sin(L*xgl(j)) + L*c*sin(L*xgl(i))*cos(L*xgl(j)))*wgl(i)*wgl(j);
        elseif prob == 2
            vecW(i+(j-1)*(N+1)) = (L^2/4)*wgl(i)*wgl(j);
        end
    end
end

W = diag(vecW);

CC = A'*W*A;

M = diag([vecW vecW]);

a_eig = sort(eig(CC,M));

noz = 0;
for i=length(a_eig):-1:1
    if abs(real(a_eig(i))) < 1e-2
        a_eig(i) = [];
        noz = noz + 1;
    end
end

a_eig = [0 ; a_eig];

no = [0 1 2 3 4 5 6 7 8 9 10];

figure(1)
plot(no,a_eig(1:11),'o',...
    'MarkerFaceColor','r',...
    'MarkerSize',5)
grid on
if prob==1
title(['Result for N=',num2str(N), ' and c=',num2str(c)])
elseif prob == 2
    title(['Result for N=',num2str(N), ' and angle=',num2str(angle)])
end
xlabel('Value of eigenvalue')
ylabel('Number of eigenvalue')
%ylim([0 10])
xlim([0 10])
                