clear all
close all
clc

WindowSize

xe(1) = 0;
ye(1) = 0;

xe(2) = 1;
ye(2) = 0;

xe(3) = 0.5;
ye(3) = sqrt(3)/2;

M = [ 1 1 1 ; xe ; ye ];

n = 10;

xi = linspace(0,1,n);
eta = linspace(0,1,n);

[Xi,Eta] = meshgrid(xi,eta);

X = 0.5+(2*Xi-1).*(1-Eta)/2;
Y = sqrt(3)/2*Eta;

for i=1:n
    for j=1:n
        Vec = [ 1 ; X(i,j) ; Y(i,j) ];
        Z = M\Vec;
        Z1(i,j) = Z(1);
        Z2(i,j) = Z(2);
        Z3(i,j) = Z(3);
    end
end

figure(win{1,3}{:})
subplot(1,3,1)
pcolor(X,Y,Z1)
shading interp
axis equal
subplot(1,3,2)
pcolor(X,Y,Z2)
shading interp
axis equal
subplot(1,3,3)
pcolor(X,Y,Z3)
shading interp
axis equal

%%


% A = det(M)/2;
% 
% MM = [ xe(2)*ye(3)-xe(3)*ye(2) ye(2)-ye(3) xe(3)-xe(2)
%        xe(3)*ye(1)-xe(1)*ye(3) ye(3)-ye(1) xe(1)-xe(3)
%        xe(1)*ye(2)-xe(2)*ye(1) ye(1)-ye(2) xe(2)-xe(1) ]/(2*A);

MM = inv(M);

%%
for i=1:n
    for j=1:n
        for E=1:3
            Ind = rem(E+1,3)+3*(E+1==3);
            Wx{E}(i,j) = MM(E,:)*[ 1 ; X(i,j) ; Y(i,j) ]*MM(Ind,2) - MM(Ind,:)*[ 1 ; X(i,j) ; Y(i,j) ]*MM(E,2);
            Wy{E}(i,j) = MM(E,:)*[ 1 ; X(i,j) ; Y(i,j) ]*MM(Ind,3) - MM(Ind,:)*[ 1 ; X(i,j) ; Y(i,j) ]*MM(E,3);
        end
    end
end


figure(win{1,3}{:})
subplot(1,3,1)
quiver(X,Y,Wx{1},Wy{1})
% shading interp
axis equal
subplot(1,3,2)
quiver(X,Y,Wx{2},Wy{2})
% shading interp
axis equal
subplot(1,3,3)
quiver(X,Y,Wx{3},Wy{3})
% shading interp
axis equal


figure(win{1,1}{:})
contour(X,Y,sqrt(Wx{1}.^2+Wy{1}.^2),'ShowText','on')
% shading interp
axis equal
colormap('jet')
colorbar



%%
for i=1:n
    for j=1:n
        for E=1:3
            Ind = rem(E+1,3)+3*(E+1==3);
            Wx{E}(i,j) = MM(E,:)*[ 1 ; X(i,j) ; Y(i,j) ]*MM(Ind,2) + MM(Ind,:)*[ 1 ; X(i,j) ; Y(i,j) ]*MM(E,2);
            Wy{E}(i,j) = MM(E,:)*[ 1 ; X(i,j) ; Y(i,j) ]*MM(Ind,3) + MM(Ind,:)*[ 1 ; X(i,j) ; Y(i,j) ]*MM(E,3);
        end
    end
end

figure(win{1,3}{:})
subplot(1,3,1)
quiver(X,Y,Wx{1},Wy{1})
% shading interp
axis equal
subplot(1,3,2)
quiver(X,Y,Wx{1},-Wy{1})
% shading interp
axis equal
subplot(1,3,3)
quiver(X,Y,-Wx{1},Wy{1})
% shading interp
axis equal