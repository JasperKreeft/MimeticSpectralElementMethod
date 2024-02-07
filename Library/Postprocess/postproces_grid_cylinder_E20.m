function [Meshp,hp,ep] = postproces_grid_cylinder_E20(R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights

global N numElements
global nn

nn = 30;%4*N; % 50;
[xip,wp] = GLLnodes(nn-1); etap = xip;
Meshp.W = kron(wp,wp)';

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xglobal = [ -1.5 -0.75 0 0.75 3 -1.5 -0.75 0 0.75 3 ...
%             -1.5 -0.75 0.75 3 ...
%             -1.5 -0.75 0 0.75 3 -1.5 -0.75 0 0.75 3 ...
%             -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 0 -sqrt(2)/4 ];
% yglobal = [ -0.75 -0.75 -0.75 -0.75 -0.75 -0.625 -0.625 -0.625 -0.625 -0.625 ...
%              0 0 0 0 ...
%              0.625 0.625 0.625 0.625 0.625 0.75 0.75 0.75 0.75 0.75 ...
%              0 -sqrt(2)/4 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 ];

xglobal = [ -2 -1 0 1 2 -2 -1 0 1 2 ...
            -2 -1 1 2 ...
            -2 -1 0 1 2 -2 -1 0 1 2 ...
            -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 0 -sqrt(2)/4 ];
yglobal = [ -2 -2 -2 -2 -2 -1 -1 -1 -1 -1 ...
             0 0 0 0 ...
             1 1 1 1 1 2 2 2 2 2 ...
             0 -sqrt(2)/4 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 ];


corners = [  1  2  7  6
             2  3  8  7
             3  4  9  8
             4  5 10  9
             6  7 12 11
             7 26 25 12
             7  8 27 26
             8  9 28 27
            28  9 13 29
             9 10 14 13
            11 12 16 15
            12 25 32 16
            32 31 17 16
            31 30 18 17
            29 13 18 30
            13 14 19 18
            15 16 21 20
            16 17 22 21
            17 18 23 22
            18 19 24 23 ];


% Global Mesh
Meshp.X      = zeros(nn^2,numElements);
Meshp.Y      = zeros(nn^2,numElements);
Meshp.dXdXi  = zeros(nn^2,numElements);
Meshp.dXdEta = zeros(nn^2,numElements);
Meshp.dYdXi  = zeros(nn^2,numElements);
Meshp.dYdEta = zeros(nn^2,numElements);
Meshp.J      = zeros(nn^2,numElements);

for i=[1:5 10 11 16:20]

[ Meshp.X(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i) ] = transfinitemapping(xglobal(corners(i,:)),nn,xip);
[ Meshp.Y(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i) ] = transfinitemapping(yglobal(corners(i,:)),nn,xip);

end


tblr = zeros(8,4);
tblr(1,4) = 1; tblr(2,1) = 1; tblr(3,1) = 1; tblr(4,3) = 1;
tblr(5,4) = 1; tblr(6,2) = 1; tblr(7,2) = 1; tblr(8,3) = 1;

theta = [ -3/4*pi     -pi
          -3/4*pi   -pi/2 
            -pi/2   -pi/4
            -pi/4       0
               pi  3/4*pi
           3/4*pi    pi/2
             pi/2    pi/4
                0    pi/4 ];

k=0;
for i=[ 6:9 12:15 ]

    k=k+1;
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];
        
[Meshp.X(:,i),Meshp.Y(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i)] = ...
                transfinitemapping_cylinder(cor,R,theta(k,:),tblr(k,:),nn,xip);
            
end

Meshp.J = Meshp.dXdXi.*Meshp.dYdEta-Meshp.dXdEta.*Meshp.dYdXi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%