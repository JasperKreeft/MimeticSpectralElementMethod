global nn

phimin = min(min(phi));
phimax = max(max(phi));

figure

subplot(1,3,1)
% surf(X,Y,PHI(globalnr_0))
dotcolorplot(PHI,Mesh.X,Mesh.Y,phimin,phimax);
colorbar off
view([0 0 1])
axis('square')
title('\phi')

subplot(1,3,2)
for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
phip = reshape(phi(:,i),nn,nn);
pcolor(Xp,Yp,phip)
hold on
end
shading interp
axis('square')
% colorbar
set(gca,'clim',[phimin phimax])
title('\phi')

subplot(1,3,3)
for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
errp = reshape(abs(phi_ex(:,i)-phi(:,i)),nn,nn);
surf(Xp,Yp,errp)
hold on
end
shading interp
axis('square')
% colorbar
title('error = |\phi_e_x-\phi_n_u_m|')