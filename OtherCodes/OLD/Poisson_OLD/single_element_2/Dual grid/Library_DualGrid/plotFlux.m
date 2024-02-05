
figure

subplot(1,3,1)
pcolor(Xp,Yp,qx)
shading interp
title('q_\xi')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

subplot(1,3,2)
pcolor(Xp,Yp,qy)
shading interp
title('q_\eta')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

subplot(1,3,3)
pcolor(Xp,Yp,qq)
shading interp
title('q')
axis('square')
% colorbar
set(gca,'clim',[0 pi])
hold on
quiver(Xp(1:2:nn-1,1:2:nn-1),Yp(1:2:nn-1,1:2:nn-1),qx(1:2:nn-1,1:2:nn-1)',qy(1:2:nn-1,1:2:nn-1)','w')