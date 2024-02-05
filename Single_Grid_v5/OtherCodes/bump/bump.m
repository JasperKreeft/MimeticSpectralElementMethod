% De bump komt van Van der Bas thesis

clear all
close all
clc

% xi = linspace(-1,1,4);
% eta = xi;




%%%%%%%%%%%%%%%%%%%%%%%%%%

rand = 'circle';

N = 10;
dx = 1/N;

% x = linspace(0,1,N+1);
x = (GLLnodes(N)+1)/2;

switch rand
    case 'naca0012'
        y = +0.6*(0.29690*sqrt(x) - 0.12600*x - 0.35160*x.^2 + 0.28430*x.^3 -0.10360*x.^4);
    case 'circle'
        R = 1; yp = sqrt(R^2-0.25); y = sqrt(R^2-(x-0.5).^2) - yp;
    otherwise
        y = zeros(1,N+1);
end

% plot(x,y)
% axis equal


%%%%%%%%%%%%%%%%%%%%%%%%%%%


X = ones(N+1,1)*x;
Y = X';

for i=0:N
    X(:,i+1) = X(:,i+1);% + y'*(N-i)/N-y'*i/N;
    Y(i+1,:) = Y(i+1,:) + y*(N-i)/N;% -y*i/N;
end

surf(X,Y,ones(N+1))
view([0 0 1])