%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_ex = sin(m*pi*Xp).*sin(m*pi*Yp);
% phi_ex = 1/4*(cos(pi*Xp)+1).*(cos(pi*Yp)+1);
% phi_ex = Xp.^2.*Yp.^2-Xp.^2-Yp.^2+1;

% surf(Xp,Yp,phi_ex)
% shading interp
% contour(Xp,Yp,phi_ex','k')
% colorbar
% set(gca,'clim',[0 1])
% set(gca,'clim',[-1 1])
% hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_ex = m*pi*cos(m*pi*Xp).*sin(m*pi*Yp);
v_ex = m*pi*sin(m*pi*Xp).*cos(m*pi*Yp);
% u_ex = -pi/4*sin(pi*Xp).*(cos(pi*Yp)+1);
% v_ex = -pi/4*(cos(pi*Xp)+1).*sin(pi*Yp);
% u_ex = 2*Xp.*Yp.^2-2*Xp;
% v_ex = 2*Xp.^2.*Yp-2*Yp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_ex = -2*m^2*pi^2*sin(m*pi*Xp).*sin(m*pi*Yp);
% f_ex = -pi^2/4*(cos(pi*Xp)*(cos(pi*Yp)+1)+(cos(pi*Xp)+1)*cos(pi*Yp));
% f_ex = 4*Xp.*Yp-2/3*Xp.^3.*Yp-2/3*Xp.*Yp.^3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%