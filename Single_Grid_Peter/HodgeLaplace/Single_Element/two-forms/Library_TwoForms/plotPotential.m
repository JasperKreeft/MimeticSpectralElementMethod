global nn

Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);

phimin = min(min(phi));
phimax = max(max(phi));

figure

subplot(1,3,1)
% % STILL EMPTY
% view([0 0 1])
% axis('square')
% title('\phi')
Phip = reshape(phi_interp,nn,nn);
surf(Xp,Yp,Phip)
shading interp
axis('square')
% colorbar
set(gca,'clim',[phimin phimax])
title('\phi_{interp}')


subplot(1,3,2)
Phip = reshape(phi,nn,nn);
surf(Xp,Yp,Phip)
shading interp
axis('square')
% colorbar
set(gca,'clim',[phimin phimax])
title('\phi')

subplot(1,3,3)
Phip = reshape(abs(phi_ex-phi),nn,nn);
surf(Xp,Yp,Phip)
shading interp
axis('square')
% colorbar
title('error = |\phi_e_x-\phi_n_u_m|')