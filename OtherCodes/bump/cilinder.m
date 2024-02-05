% Idee van Kwakkel thesis fig. 5.18

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cilinder

Rc = 1;
xc = linspace(-1,1,100); % coordinates of cilinder
yc = sqrt(Rc^2-xc.^2);

plot(xc,yc,'r')
axis equal
axis([-4 4 0 2])
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh around cilinder

theta = 180:-30:0;
dR = 0.5*(2/sind(60)-1);

xv =zeros(3,length(theta)); yv = xv;
for i=1:3
Rv = Rc+(i-1)*dR;
xv(i,:) = Rv*cosd(theta); % coordinates of vertices of mesh
yv(i,:) = Rv*sind(theta);
end
yv = min(yv,2);

hold on
plot(xv,yv,'bx')
plot(xv(3,:),yv(3,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outer mesh

xol = [-4 -3];
yol = [0 (Rc+2*dR)*sind(30) 2];
[Xol,Yol] = meshgrid(xol,yol);
Xol(3,2)= -2.5;
plot(Xol,Yol,'gx')

xor = -xol;
yor = yol;
[Xor,Yor] = meshgrid(xor,yor);
Xor(3,2)= 2.5;
plot(Xor,Yor,'gx')

% close