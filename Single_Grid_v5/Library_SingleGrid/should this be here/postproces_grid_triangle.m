function [Meshp,hp,ep] = postproces_grid_triangle()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 

global N nn xi

nn = 8*N; % 50;N+1;%
[xip,wp] = GLLnodes(nn-1);
Meshp.W = wp'*wp;

[hp ep] = MimeticpolyVal(xip,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = sqrt(3)/3;

xglobal = [ -1  0  1 1/2 0  -1/2 0 ];
yglobal = [ -a -a -a a/2 2*a a/2 0 ];

corners = [  1 2 7 6
             2 3 4 7
             7 4 5 6 ];


% Global Mesh
Meshp.X      = zeros(nn^2,3);
Meshp.Y      = zeros(nn^2,3);
Meshp.J      = zeros(nn^2,3);
Meshp.dXdXi  = zeros(nn^2,3);
Meshp.dXdEta = zeros(nn^2,3);
Meshp.dYdXi  = zeros(nn^2,3);
Meshp.dYdEta = zeros(nn^2,3);

for i=1:3

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