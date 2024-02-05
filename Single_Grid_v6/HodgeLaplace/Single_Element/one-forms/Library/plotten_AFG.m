global nn

Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);

XYlim = [min(min(Xp)) max(max(Xp)) min(min(Yp)) max(max(Yp))];

%%

sp = reshape(ssigma,nn,nn);
up = reshape(uu,nn,nn);
vp = reshape(vv,nn,nn);
ng_sxp = reshape(ng_sx,nn,nn);
ng_syp = reshape(ng_sy,nn,nn);
div_uvp = reshape(div_uv,nn,nn);
if exist('ffx','var')
fxp = reshape(ffx,nn,nn);
fyp = reshape(ffy,nn,nn);
end

%%

figure
subplot(2,3,1)
pcolor(Xp,Yp,up)
shading interp
% colorbar
title('u')
axis equal
axis(XYlim)
subplot(2,3,2)
pcolor(Xp,Yp,vp)
shading interp
% colorbar
title('v')
axis equal
axis(XYlim)
subplot(2,3,3)
pcolor(Xp,Yp,div_uvp)
shading interp
% colorbar
title('div u')
axis equal
axis(XYlim)
subplot(2,3,4)
pcolor(Xp,Yp,sp)
shading interp
% set(gca,'clim',[-5 5])
% colorbar
title('\sigma')
axis equal
axis(XYlim)
subplot(2,3,5)
pcolor(Xp,Yp,ng_sxp)
% set(gca,'clim',[-5 5])
shading interp
% colorbar
title('curl \sigma (x)')
axis equal
axis(XYlim)
subplot(2,3,6)
pcolor(Xp,Yp,ng_syp)
shading interp
% colorbar
title('curl \sigma (y)')
axis equal
axis(XYlim)

% %%
figure
% subplot(2,3,1)
% contourf(Xp,Yp,up,20)
% colorbar
% title('u')
% axis equal
% axis(XYlim)
% subplot(2,3,2)
% contourf(Xp,Yp,vp,20)
% colorbar
% title('v')
% axis equal
% axis(XYlim)
% subplot(2,3,3)
% contourf(Xp,Yp,wp,-5:5)
% colorbar
% title('w')
% axis equal
% axis(XYlim)
% subplot(2,3,4)
% contourf(Xp,Yp,Pp,-5:5)
% colorbar
% title('p')
% axis equal
% axis(XYlim)
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
end
