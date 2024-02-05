phimin = -1; % min(phi);
phimax =  1; % max(phi);

figure

subplot(1,3,1)
% surf(XEGEG,YEGEG,PHI)
% surf(XGG,YGG,PHI)
% dotcolorplot(PHI,XEGEG,YEGEG,phimin,phimax)
dotcolorplot(PHI,XGG,YGG,phimin,phimax)
colorbar off
view([0 0 1])
axis('square')
title('\phi')


subplot(1,3,2)
pcolor(Xp,Yp,phi)
shading interp
axis('square')
% colorbar
set(gca,'clim',[phimin phimax])
title('\phi')


subplot(1,3,3)
surf(Xp,Yp,abs(phi_ex-phi))
shading interp
axis('square')
% colorbar
% set(gca,'clim',[phimin phimax])
title('error = |\phi_e_x-\phi_n_u_m|')