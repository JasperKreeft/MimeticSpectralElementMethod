function [Meshp,hp,ep] = postproces_grid_naca0012_E16(R)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights

global N numElements
global nn

% close all
% R = 0.5; N=10;

nn = N+1;
[xip,wp] = GLLnodes(nn-1); etap = xip;
Meshp.W = kron(wp,wp)';

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


xglobal = [ -2 0 2 2 2 0 -2 -2 ...
            -1 0 1 1 1 0 -1 -1 ...
            -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 0 -sqrt(2)/4 -0.5 ];
yglobal = [ -2 -2 -2 0 2 2 2 0 ...
            -1 -1 -1 0 1 1 1 0 ...
            -sqrt(2)/4 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4  0 ];
         

% plot(xglobal,yglobal,'o','markerface','b')
hold on
axis equal
axis([-2 2 -2 2])


corners = [  1  2 10  9 % outer-ring
             2  3 11 10
             3  4 12 11
             4  5 13 12
             5  6 14 13
             6  7 15 14
             7  8 16 15
             8  1  9 16
             9 10 18 17 % inner-ring
            10 11 19 18
            11 12 20 19
            12 13 21 20
            13 14 22 21
            14 15 23 22
            15 16 24 23
            16  9 17 24 ];


% Global Mesh
Meshp.X      = zeros(nn^2,16);
Meshp.Y      = zeros(nn^2,16);
Meshp.dXdXi  = zeros(nn^2,16);
Meshp.dXdEta = zeros(nn^2,16);
Meshp.dYdXi  = zeros(nn^2,16);
Meshp.dYdEta = zeros(nn^2,16);
Meshp.J      = zeros(nn^2,16);

for i=1:8

[ Meshp.X(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i) ] = transfinitemapping(xglobal(corners(i,:)),nn,xip);
[ Meshp.Y(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i) ] = transfinitemapping(yglobal(corners(i,:)),nn,xip);

% pcolor(reshape(Meshp.X(:,i),nn,nn),reshape(Meshp.Y(:,i),nn,nn),(i-1)/15*ones(nn))

end


tblr = zeros(8,4);
tblr(1,1) = 1; tblr(2,1) = 1; tblr(3,1) = 1; tblr(4,1) = 1;
tblr(5,1) = 1; tblr(6,1) = 1; tblr(7,1) = 1; tblr(8,1) = 1;

theta = [ -3/4*pi   -pi/2 
            -pi/2   -pi/4
            -pi/4       0
                0  1/4*pi
           1/4*pi    pi/2
             pi/2    3*pi/4
           3*pi/4    pi
              -pi   -3/4*pi ];

k=0;
for i=9:16

    k=k+1;
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];

[Meshp.X(:,i),Meshp.Y(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i)] = ...
                transfinitemapping_naca0012(cor,R,theta(k,:),tblr(k,:),nn,xip);            

% pcolor(reshape(Meshp.X(:,i),nn,nn),reshape(Meshp.Y(:,i),nn,nn),(i-1)/15*ones(nn))

end

Meshp.J = Meshp.dXdXi.*Meshp.dYdEta-Meshp.dXdEta.*Meshp.dYdXi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%