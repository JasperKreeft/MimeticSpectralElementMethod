clc;
clear all;
close all;

% This program solves the curl-curl eigenvalue problem on a rectangular
% domain using nodal spectral elements. The purpose of this code is to
% compare the results with the mimetic spectral element approach.
%
% 19 August 2010, Marc Gerritsma

P_range = 2%:2:30;

for P=P_range
P
% P=10;        % Polynomial degree
L=pi;        % Size of the domain
c = 0.2;     % Curvature of the gridlines


[xgl wgl] = lglnodes(P) ;  % GL-points and weights
xgl = flipud(xgl)'      ;
wgl = flipud(wgl)'      ;

%% Legendre polynomial of degree P evaluated at the Gauss-Lobatto points

for i=1:P+1
    Lprev = 1;
    Lthis = xgl(i);
    for k=1:P-1;
        LN(i) = ( (2*k+1)*xgl(i)*Lthis - k*Lprev )/(k+1);
        Lprev = Lthis;
        Lthis = LN(i);
    end
end

% Setting up the derivative matrix

for i=1:P+1
    for j=1:P+1
        if i==j
            if i==1
                dh(i,j) = -P*(P+1)/4;
            elseif i==P+1
                dh(i,j) = P*(P+1)/4;
            else
                dh(i,j) = 0;
            end
        else
            dh(i,j) = LN(j)/( LN(i)*(xgl(j)-xgl(i)));
        end
    end
end

%  Loop over the basis functions
A = zeros((P+1)*(P+1),2*(P+1)*(P+1));
for i=1:P+1
    for j=1:P+1
        % Loop over the integration points
        for p=1:P+1
            for q=1:P+1
                dxdxi  = (L/2)*(1 + L*c*cos(L*xgl(p))*sin(L*xgl(q)));
                dxdeta = (L/2)*L*c*sin(L*xgl(p))*cos(L*xgl(q));
                dydxi  = (L/2)*L*c*cos(L*xgl(p))*sin(L*xgl(q));
                dydeta = (L/2)*(1 + L*c*sin(L*xgl(p))*cos(L*xgl(q)));
                J = dxdxi*dydeta - dxdeta*dydxi;
                eqnr = p + (q-1)*(P+1);
                if q==j
                    var  = i + (q-1)*(P+1);
                    A(eqnr,var) = A(eqnr,var) + dydeta*dh(i,p)/J;  % dv/dx :  dv/dxi*dxi/dx = dx/dxi * dy/deta/J
                    var = (P+1)*(P+1) + i + (q-1)*(P+1);
                    A(eqnr,var) = A(eqnr,var) + dxdeta*dh(i,p)/J;
                end
                if p==i
                    var = (P+1)*(P+1) + p + (j-1)*(P+1);
                    A(eqnr,var) = A(eqnr,var) - dxdxi*dh(j,q)/J;
                    var = p + (j-1)*(P+1);
                    A(eqnr,var) = A(eqnr,var) - dydxi*dh(j,q)/J;
                end
            end
        end
    end
end

for i=1:P+1
    for j=1:P+1
        vecW(i+(j-1)*(P+1)) = (L^2/4)*(1 + L*c*cos(L*xgl(i))*sin(L*xgl(j)) + L*c*sin(L*xgl(i))*cos(L*xgl(j)))*wgl(i)*wgl(j);
    end
end

W = diag(vecW);

CC = A'*W*A;

M = diag([vecW vecW]);

a_eig = sort(eig(CC,M));

noz = 0;
for i=length(a_eig):-1:1
    if abs(real(a_eig(i))) < 1e-8
        a_eig(i) = [];
        noz = noz + 1;
    end
end

% keyboard

ind = min(length(a_eig),10);
a_eig_total(1:ind,P) = a_eig(1:ind);
exact = [1 1 2 4 4 5 5 8 9 9]';
error(1:ind,P) = abs(a_eig_total(1:ind,P)-exact(1:ind));
clearallbut a_eig_total error P_range c




% no = [0 1 2 3 4 5 6 7 8 9 10];
% 
% figure(1)
% plot(no(1:ind+1),[0 ; a_eig(1:ind)],'o',...
%     'MarkerFaceColor','r',...
%     'MarkerSize',5)
% grid on
% title(['Result for N=',num2str(P)])
% xlabel('Value of eigenvalue')
% ylabel('Number of eigenvalue')
% axis([0 10 0 10])
% 
% pause
end

figure(2)
semilogy(P_range,error(1,P_range),'-ob','markerface','b')
hold on
semilogy(P_range,error(2,P_range),'-or','markerface','r')
semilogy(P_range,error(3,P_range)','-og','markerface','g')
semilogy(P_range,error(4,P_range)','-ok','markerface','k')
semilogy(P_range,error(5,P_range),'-om','markerface','m')
semilogy(P_range,error(6,P_range)','-oc','markerface','c')
semilogy(P_range,error(7,P_range),'-oy','markerface','y')
semilogy(P_range,error(8,P_range)','-^b','markerface','b')
semilogy(P_range,error(9,P_range)','-^r','markerface','r')
grid on
legend('1','2','3','4','5','6','7','8','9',1)
% axis([0 20 1e-10 1e2])
xlabel('P')
ylabel('error eigenvalues')
title(['Convergence of first ten eigenvalues for c=' num2str(c)])