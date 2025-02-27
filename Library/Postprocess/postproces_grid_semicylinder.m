function [Meshp,hp,ep] = postproces_grid_semicylinder(R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights

global N numElements
global nn

nn = 8*N; % 50;
[xip,wp] = GLLnodes(nn-1); etap = xip;
Meshp.W = kron(wp,wp)';

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xglobal = [ -1.5 -0.75 0 0.75 3 -1.5 -0.75 0.75 3 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 ];
yglobal = [ -0.75 -0.75 -0.75 -0.75 -0.75 0 0 0 0 0 -sqrt(2)/4 -0.5 -sqrt(2)/4 0 ];


corners = [  1  2  7  6
             2 11 10  7
             2  3 12 11
             3  4 13 12
            13  4  8 14
             4  5  9  8 ];


% Global Mesh
Meshp.X      = zeros(nn^2,numElements);
Meshp.Y      = zeros(nn^2,numElements);
Meshp.dXdXi  = zeros(nn^2,numElements);
Meshp.dXdEta = zeros(nn^2,numElements);
Meshp.dYdXi  = zeros(nn^2,numElements);
Meshp.dYdEta = zeros(nn^2,numElements);
Meshp.J      = zeros(nn^2,numElements);

for i=[1 6]

    [ Meshp.X(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i) ] = transfinitemapping(xglobal(corners(i,:)),nn,xip);
    [ Meshp.Y(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i) ] = transfinitemapping(yglobal(corners(i,:)),nn,xip);

end


tblr = zeros(4,4);
tblr(1,4) = 1; tblr(2,1) = 1; tblr(3,1) = 1; tblr(4,3) = 1;

theta = [ -3/4*pi     -pi
          -3/4*pi   -pi/2 
            -pi/2   -pi/4
            -pi/4       0 ];

k=0;
for i=2:5

    k=k+1;
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];

[Meshp.X(:,i),Meshp.Y(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i)] = ...
                transfinitemapping_cylinder(cor,R,theta(k,:),tblr(k,:),nn,xip);

end

Meshp.J = Meshp.dXdXi.*Meshp.dYdEta-Meshp.dXdEta.*Meshp.dYdXi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%