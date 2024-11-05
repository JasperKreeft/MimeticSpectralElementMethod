function plotFlux(Meshp,qx,qy,qMag,qx_int,qy_int,q_interp,FunctionType,Domain)

global nn numElements

figure

Xplim = [min(min(Meshp.X)) max(max(Meshp.X))];
Yplim = [min(min(Meshp.Y)) max(max(Meshp.Y))];

for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
Qpx = reshape(qx(:,i),nn,nn);
subplot(2,3,1)
pcolor(Xp,Yp,Qpx)
shading interp
title('q_x')
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[-pi pi])
end

subplot(2,3,2)
for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
Qpy = reshape(qy(:,i),nn,nn);
pcolor(Xp,Yp,Qpy)
hold on
end
shading interp
title('q_y')
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[-pi pi])

subplot(2,3,3)
for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
Qp = reshape(qMag(:,i),nn,nn);
pcolor(Xp,Yp,Qp)
hold on
end
shading interp
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[0 pi])
title('absolute flux')
hold on
% quiver(Xp(1:2:nn-1,1:2:nn-1),Yp(1:2:nn-1,1:2:nn-1),Qpx(1:2:nn-1,1:2:nn-1),Qpy(1:2:nn-1,1:2:nn-1),'w')

if exist('qx_int','var')
subplot(2,3,4)
for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
Qpx_int = reshape(qx_int(:,i),nn,nn);
pcolor(Xp,Yp,Qpx_int)
hold on
end
shading interp
title('q_x interp')
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[-pi pi])

subplot(2,3,5)
for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
Qpy_int = reshape(qy_int(:,i),nn,nn);
pcolor(Xp,Yp,Qpy_int)
hold on
end
shading interp
title('q_y interp')
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[-pi pi])

subplot(2,3,6)
for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
Qp = reshape(q_interp(:,i),nn,nn);
pcolor(Xp,Yp,Qp)
hold on
end
shading interp
axis equal
axis([Xplim Yplim])
% colorbar
% set(gca,'clim',[0 pi])
title('absolute flux interp')
hold on
quiver(Xp(1:2:nn-1,1:2:nn-1),Yp(1:2:nn-1,1:2:nn-1),Qpx_int(1:2:nn-1,1:2:nn-1),Qpy_int(1:2:nn-1,1:2:nn-1),'w')
end



if strcmp(FunctionType,'nozzle') && strcmp(Domain,'LavalNozzle')
    figure
    for i=1:numElements
    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
    Qp = reshape(qMag(:,i),nn,nn);
    pcolor(Xp,Yp,Qp)
    hold on
    end
    shading interp
    axis equal
    axis([Xplim Yplim])
    title('absolute flux')
    hold on
    quiver(Xp(1:2:nn-1,1:2:nn-1),Yp(1:2:nn-1,1:2:nn-1),Qpx(1:2:nn-1,1:2:nn-1),Qpy(1:2:nn-1,1:2:nn-1),'w')
    colorbar
end