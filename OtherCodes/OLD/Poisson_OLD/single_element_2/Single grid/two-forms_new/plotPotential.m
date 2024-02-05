phimin = min(min(phi));
phimax = max(max(phi));

figure

subplot(1,3,1)
% STILL EMPTY
view([0 0 1])
axis('square')
title('\phi')

subplot(1,3,2)
surf(Xp,Yp,phi)
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
title('error = |\phi_e_x-\phi_n_u_m|')