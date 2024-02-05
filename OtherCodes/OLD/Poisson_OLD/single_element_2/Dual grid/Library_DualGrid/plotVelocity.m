umin    = -pi;
umax    =  pi;
velomin =   0;
velomax =  pi;

figure

subplot(1,3,1)
pcolor(Xp,Yp,ux)
shading interp
title('u')
axis('square')
% colorbar
set(gca,'clim',[umin umax])

subplot(1,3,2)
pcolor(Xp,Yp,uy)
shading interp
title('v')
axis('square')
% colorbar
set(gca,'clim',[umin umax])

subplot(1,3,3)
pcolor(Xp,Yp,velo)
shading interp
axis('square')
% colorbar
set(gca,'clim',[velomin velomax])
title('absolute velocity')
hold on
quiver(Xp(1:2:nn-1,1:2:nn-1),Yp(1:2:nn-1,1:2:nn-1),ux(1:2:nn-1,1:2:nn-1)',uy(1:2:nn-1,1:2:nn-1)','w')
