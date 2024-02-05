figure

Xplim = [min(min(Xp)) max(max(Xp))];
Yplim = [min(min(Yp)) max(max(Yp))];

subplot(2,3,1)
pcolor(Xp,Yp,qx)
shading interp
title('q_x')
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[-pi pi])

subplot(2,3,2)
pcolor(Xp,Yp,qy)
shading interp
title('q_y')
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[-pi pi])

subplot(2,3,3)
pcolor(Xp,Yp,velo)
shading interp
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[0 pi])
title('absolute flux')
hold on
quiver(Xp(1:2:nn-1,1:2:nn-1),Yp(1:2:nn-1,1:2:nn-1),qx(1:2:nn-1,1:2:nn-1),qy(1:2:nn-1,1:2:nn-1),'w')

if exist('qx_int','var')
subplot(2,3,4)
pcolor(Xp,Yp,qx_int)
shading interp
title('q_x interp')
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[-pi pi])

subplot(2,3,5)
pcolor(Xp,Yp,qy_int)
shading interp
title('q_y interp')
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[-pi pi])

subplot(2,3,6)
pcolor(Xp,Yp,q_interp)
shading interp
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[0 pi])
title('absolute flux interp')
hold on
quiver(Xp(1:2:nn-1,1:2:nn-1),Yp(1:2:nn-1,1:2:nn-1),qx_int(1:2:nn-1,1:2:nn-1),qy_int(1:2:nn-1,1:2:nn-1),'w')
end



if strcmp(FunctionType,'nozzle') && strcmp(Domain,'LavalNozzle')
    figure
    pcolor(Xp,Yp,velo)
    shading interp
    axis equal
    axis([Xplim Yplim])
    title('absolute flux')
    hold on
    quiver(Xp(1:2:nn-1,1:2:nn-1),Yp(1:2:nn-1,1:2:nn-1),qx(1:2:nn-1,1:2:nn-1),qy(1:2:nn-1,1:2:nn-1),'w')
    colorbar
end