function [Meshp,hp,ep] = postproces_grid_annulus(R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights

global N numElements
global nn

nn = N+1; % 50;
[xip,wp] = GLLnodes(nn-1);
xip = linspace(-1,1,nn);
etap = xip;
Meshp.W = kron(wp,wp)';

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xglobal = [ 0 1 0 -1 0 R 0 -R ];
% yglobal = [ -1 0 1 0 -R 0 R 0 ];

a = 1/2*sqrt(2);
xglobal = [ a a -a -a R*a R*a -R*a -R*a ];
yglobal = [ -a a a -a -R*a R*a R*a -R*a ];

corners = [ 1 2 6 5
            2 3 7 6
            3 4 8 7
            4 1 5 8 ];


% Global Mesh
Meshp.X      = zeros(nn^2,numElements);
Meshp.Y      = zeros(nn^2,numElements);
Meshp.dXdXi  = zeros(nn^2,numElements);
Meshp.dXdEta = zeros(nn^2,numElements);
Meshp.dYdXi  = zeros(nn^2,numElements);
Meshp.dYdEta = zeros(nn^2,numElements);
Meshp.J      = zeros(nn^2,numElements);

% theta = [ -pi/2      0   
%               0   pi/2
%            pi/2     pi
%             -pi -1/2*pi ];
        
theta = [ -pi/4    pi/4   
           pi/4    3*pi/4
         3*pi/4    5*pi/4
         -3*pi/4   -pi/4 ];

for i=1:4

    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];
        
[Meshp.X(:,i),Meshp.Y(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i)] = ...
                transfinitemapping_annulus_post(cor,[1 R],theta(i,:),nn,xip);
end

Meshp.J = Meshp.dXdXi.*Meshp.dYdEta-Meshp.dXdEta.*Meshp.dYdXi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%