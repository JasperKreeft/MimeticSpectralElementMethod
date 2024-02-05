global c m

[Dp,Gd] = topology(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 
nn = 100;
[xi,wg] = Gnodes(nn); eta = xi;
Xi = xi'*ones(1,nn);  Eta = Xi';
Wgg = wg'*wg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hEG dhEGdx] = LagrangeVal(xi,N,3);
hG = LagrangeVal(xi,N,2);
eEG = EdgeVal(dhEGdx);

switch gridtype
    case 'standard'
        Xic = Xi;
        Etac = Eta;
    case 'sinecurve'
        Xic = Xi+c*sin(pi*Xi).*sin(pi*Eta);
        Etac = Eta+c*sin(pi*Xi).*sin(pi*Eta);
end
        
dxdxi  = 1+pi*c*cos(pi*Xi).*sin(pi*Eta);
dxdeta = pi*c*sin(pi*Xi).*cos(pi*Eta);
dydxi  = pi*c*cos(pi*Xi).*sin(pi*Eta);
dydeta = 1+pi*c*sin(pi*Xi).*cos(pi*Eta);

J = dxdxi.*dydeta-dxdeta.*dydxi;

%% Exact %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problem
    case 'sine'
        phi_ex    = sin(m*pi*Xic).*sin(m*pi*Etac);
    case 'cosine'
        phi_ex    = (0.5*cos(pi*Xic)+0.5).*(0.5*cos(pi*Etac)+0.5);
end


%% Add boundary conditions to solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = [zeros(N+2,1) [zeros(1,N); reshape(phi_in,N,N); zeros(1,N)] zeros(N+2,1)];

PHI = reshape(phi,[],1);

figure(9)
subplot(2,3,1)
surf(XEGEG,YEGEG,phi)
% view([0 0 1])
axis('square')
title('\phi')

grad_phi = Gd*phi_in;
U = grad_phi(1:N*(N+1),1);
V = grad_phi(N*(N+1)+1:2*N*(N+1),1);

u = [zeros(N+1,1) reshape(U,N+1,N) zeros(N+1,1)];
v = [zeros(1,N+1); reshape(V,N,N+1); zeros(1,N+1)];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pphi = hEG'*phi*hEG;

figure(9)
subplot(2,3,2)
pcolor(Xic,Etac,pphi)
shading interp
axis('square')
% colorbar
set(gca,'clim',[-1 1])
title('\phi')

figure(9)
subplot(2,3,3)
surf(Xic,Etac,abs(phi_ex-pphi))
shading interp
axis('square')
% colorbar
% set(gca,'clim',[-1 1])
title('error = |\phi_e_x-\phi_n_u_m|')

%% Gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uxi      = eEG'*u*hEG;
ueta     = hEG'*v*eEG;
ux = (dydeta.*uxi-dydxi.*ueta)./J;
uy = (dxdxi.*ueta-dxdeta.*uxi)./J;

figure(9)
subplot(2,3,4)
pcolor(Xic,Etac,ux)
shading interp
title('u')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

figure(9)
subplot(2,3,5)
pcolor(Xic,Etac,uy)
shading interp
title('v')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

velo = sqrt(ux.^2+uy.^2);

figure(9)
subplot(2,3,6)
pcolor(Xic,Etac,velo)
shading interp
axis('square')
% colorbar
set(gca,'clim',[0 pi])
title('absolute velocity')
hold on
quiver(Xic(1:2:nn-1,1:2:nn-1),Etac(1:2:nn-1,1:2:nn-1),ux(1:2:nn-1,1:2:nn-1)',uy(1:2:nn-1,1:2:nn-1)','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%