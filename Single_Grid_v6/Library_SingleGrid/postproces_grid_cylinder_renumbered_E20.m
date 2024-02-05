function [Meshp,hp,ep] = postproces_grid_cylinder_renumbered_E20(R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights

global N numElements
global nn

nn = 2*N; % 50;
[xip,wp] = GLLnodes(nn-1); etap = xip;
Meshp.W = kron(wp,wp)';

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xglobal = [ -1.5 -0.75 0 0.75 3 -1.5 -0.75 0.75 3 -1.5 -0.75 0 0.75 3 ...
            -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 0 -sqrt(2)/4 ...
            -1.5 -0.75 0 0.75 3 -1.5 -0.75 0 0.75 3 ];
        
yglobal = [ -0.60 -0.60 -0.60 -0.60 -0.60 0 0 0 0 0.60 0.60 0.60 0.60 0.60 ...
            0 -sqrt(2)/4 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 ...
            -0.75 -0.75 -0.75 -0.75 -0.75 0.75 0.75 0.75 0.75 0.75 ];


corners = [  6  1  2  7
            10  6  7 11
             7  2 16 15
            11  7 15 22
             2  3 17 16
             3  4 18 17
            12 11 22 21
            13 12 21 20
             4  8 19 18
             8 13 20 19
             5  9  8  4
             9 14 13  8
             1 23 24  2  % 23 24  2  1
            24 25  3  2
            25 26  4  3
            26 27  5  4
            29 28 10 11
            30 29 11 12
            31 30 12 13
            32 31 13 14 ];


% Global Mesh
Meshp.X      = zeros(nn^2,numElements);
Meshp.Y      = zeros(nn^2,numElements);
Meshp.dXdXi  = zeros(nn^2,numElements);
Meshp.dXdEta = zeros(nn^2,numElements);
Meshp.dYdXi  = zeros(nn^2,numElements);
Meshp.dYdEta = zeros(nn^2,numElements);
Meshp.J      = zeros(nn^2,numElements);

for i=[1 2 11:20]

[ Meshp.X(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i) ] = transfinitemapping(xglobal(corners(i,:)),nn,xip);
[ Meshp.Y(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i) ] = transfinitemapping(yglobal(corners(i,:)),nn,xip);

end


tblr = zeros(8,4);
tblr(:,1) = 1;

theta = [     -pi -3/4*pi
          -3/4*pi   -pi/2 
            -pi/2   -pi/4
            -pi/4       0
                0    pi/4
             pi/4    pi/2
             pi/2  3/4*pi
           3/4*pi      pi ];

R = 0.5;

k=0;
for i=[ 3 5 6 9 10 8 7 4 ]

    k=k+1;
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];
        
[Meshp.X(:,i),Meshp.Y(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i)] = ...
                transfinitemapping_cylinder(cor,R,theta(k,:),tblr(k,:),nn,xip);
            
end

Meshp.J = Meshp.dXdXi.*Meshp.dYdEta-Meshp.dXdEta.*Meshp.dYdXi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%