function [Meshp,hp,ep] = postproces_grid_triangle_E9()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 

global N nn xi

nn = 8*N; % 50;N+1;%
[xip,wp] = GLLnodes(nn-1);
Meshp.W = wp'*wp;

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = sqrt(3)/3;
b = 1/10;

xglobal = [ -1 -1+b  -1+b   -1+b/2       0 0 1-b  1-b     1 1-b/2       -1/2 1/2   -b/2          0       b/2         0 ];
yglobal = [ -a  -a  a*(b-1) a*(3/2*b-1) -a 0 -a  a*(b-1) -a a*(3/2*b-1)  a/2 a/2 a*(2-3/2*b) 2*a*(1-b) a*(2-3/2*b) 2*a ];

corners = [  1  2  3  4
             2  5  6  3
             5  7  8  6
             7  9 10  8
             3  6 11  4
             8 10 12  6
            11  6 14 13
             6 12 15 14
            14 15 16 13 ];


% Global Mesh
Meshp.X      = zeros(nn^2,9);
Meshp.Y      = zeros(nn^2,9);
Meshp.J      = zeros(nn^2,9);
Meshp.dXdXi  = zeros(nn^2,9);
Meshp.dXdEta = zeros(nn^2,9);
Meshp.dYdXi  = zeros(nn^2,9);
Meshp.dYdEta = zeros(nn^2,9);

for i=1:9

[Meshp.X(:,i),Meshp.dXdXi(:,i),Meshp.dXdEta(:,i)] = ...
                              transfinitemapping(xglobal(corners(i,:)),nn,xip);
[Meshp.Y(:,i),Meshp.dYdXi(:,i),Meshp.dYdEta(:,i)] = ...
                              transfinitemapping(yglobal(corners(i,:)),nn,xip);



% pcolor(reshape(Meshp.X(:,i),nn,nn),reshape(Meshp.Y(:,i),nn,nn),(i-1)/2*ones(nn))
% hold on
% axis equal
% axis([-1.2 1.2 -.6 1.2])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Meshp.J(:,i) = Meshp.dXdXi(:,i).*Meshp.dYdEta(:,i)-Meshp.dXdEta(:,i).*Meshp.dYdXi(:,i);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%