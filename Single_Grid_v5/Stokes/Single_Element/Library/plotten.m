global nn

Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);

XYlim = [min(min(Xp)) max(max(Xp)) min(min(Yp)) max(max(Yp))];

%%

wp = reshape(ww,nn,nn);
up = reshape(uu,nn,nn);
vp = reshape(vv,nn,nn);
Vp = reshape(velo,nn,nn);
Pp = reshape(pp,nn,nn);
if exist('ffx','var')
fxp = reshape(ffx,nn,nn);
fyp = reshape(ffy,nn,nn);
end
if exist('gg','var')
    gp = reshape(gg,nn,nn);
end

%%

% figure
% subplot(2,3,1)
% pcolor(Xp,Yp,up)
% shading interp
% % colorbar
% title('u')
% axis equal
% axis(XYlim)
% subplot(2,3,2)
% pcolor(Xp,Yp,vp)
% shading interp
% % colorbar
% title('v')
% axis equal
% axis(XYlim)
% subplot(2,3,3)
% pcolor(Xp,Yp,Vp)
% shading interp
% % colorbar
% title('abs(velo)')
% axis equal
% axis(XYlim)
% subplot(2,3,4)
% pcolor(Xp,Yp,wp)
% shading interp
% % set(gca,'clim',[-5 5])
% % colorbar
% title('w')
% axis equal
% axis(XYlim)
% subplot(2,3,5)
% pcolor(Xp,Yp,Pp)
% % set(gca,'clim',[-5 5])
% shading interp
% % colorbar
% title('p')
% axis equal
% axis(XYlim)
% subplot(2,3,6)
% % contourf(Xp,Yp,??,20)
% % colorbar
% title('??')
% axis equal
% axis(XYlim)
% 
% %%
% 
% figure
% pcolor(Xp,Yp,Vp)
% shading interp
% hold on
% nn = size(Xp,1);
% quiver(Xp(1:4:nn-1,1:4:nn-1),Yp(1:4:nn-1,1:4:nn-1),up(1:4:nn-1,1:4:nn-1),vp(1:4:nn-1,1:4:nn-1),'w')
% title('Lid-Driven Cavity Stokes')
% axis equal
% axis(XYlim)

%%
figure
subplot(2,3,1)
contourf(Xp,Yp,up,20)
colorbar
title('u')
axis equal
axis(XYlim)
subplot(2,3,2)
contourf(Xp,Yp,vp,20)
colorbar
title('v')
axis equal
axis(XYlim)
subplot(2,3,3)
contourf(Xp,Yp,wp,20)
colorbar
title('w')
axis equal
axis(XYlim)
subplot(2,3,4)
contourf(Xp,Yp,Pp,20)
colorbar
title('p')
axis equal
axis(XYlim)
if exist('ffx','var')
subplot(2,3,5)
contourf(Xp,Yp,fxp,20)
colorbar
title('f_x')
axis equal
axis(XYlim)
subplot(2,3,6)
contourf(Xp,Yp,fyp,20)
colorbar
title('f_y')
axis equal
axis(XYlim)

if exist('gg','var')
    figure
    contourf(Xp,Yp,gp,20)
    colorbar
    title('g')
    axis equal
    axis(XYlim)
end

end
