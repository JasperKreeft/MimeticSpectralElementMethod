global c m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 
nn = 100;
[xi,wg] = Gnodes(nn); eta = xi;
Xi = xi'*ones(1,nn);  Eta = Xi';
Wgg = wg'*wg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hEG eEG] = MimeticpolyVal(xi,N,3);
hG = LagrangeVal(xi,N,2);


switch gridtype
    case 'standard'
        X = Xi;
        Y = Eta;
    case 'sinecurve'
        X = Xi+c*sin(pi*Xi).*sin(pi*Eta);
        Y = Eta+c*sin(pi*Xi).*sin(pi*Eta);
end
        
dxdxi  = 1+pi*c*cos(pi*Xi).*sin(pi*Eta);
dxdeta = pi*c*sin(pi*Xi).*cos(pi*Eta);
dydxi  = pi*c*cos(pi*Xi).*sin(pi*Eta);
dydeta = 1+pi*c*sin(pi*Xi).*cos(pi*Eta);

J = dxdxi.*dydeta-dxdeta.*dydxi;

%% Exact %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch problem
    case 'sine'
        phi_ex    = sin(m*pi*X).*sin(m*pi*Y);
    case 'cosine'
        phi_ex    = (0.5*cos(pi*X)+0.5).*(0.5*cos(pi*Y)+0.5);
end

%% Add boundary conditions to solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = [zeros(N+2,1) [zeros(1,N); reshape(phi_in,N,N); zeros(1,N)] zeros(N+2,1)];

PHI = reshape(phi,[],1);

figure(9)
subplot(2,3,1)
surf(XEGEG,YEGEG,phi)
view([0 0 1])
axis('square')
title('\phi')

grad_phi = Gd*phi_in;
U = grad_phi(1:N*(N+1),1);
V = grad_phi(N*(N+1)+1:2*N*(N+1),1);

u = [zeros(N+1,1) reshape(U,N+1,N) zeros(N+1,1)];
v = [zeros(1,N+1); reshape(V,N+1,N)'; zeros(1,N+1)];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pphi = hEG'*phi*hEG;

figure(9)
subplot(2,3,2)
pcolor(X,Y,pphi)
shading interp
axis('square')
% colorbar
set(gca,'clim',[-1 1])
title('\phi')

figure(9)
subplot(2,3,3)
surf(X,Y,abs(phi_ex-pphi))
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
pcolor(X,Y,ux)
shading interp
title('u')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

figure(9)
subplot(2,3,5)
pcolor(X,Y,uy)
shading interp
title('v')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

velo = sqrt(ux.^2+uy.^2);

figure(9)
subplot(2,3,6)
pcolor(X,Y,velo)
shading interp
axis('square')
% colorbar
set(gca,'clim',[0 pi])
title('absolute velocity')
hold on
quiver(X(1:2:nn-1,1:2:nn-1),Y(1:2:nn-1,1:2:nn-1),ux(1:2:nn-1,1:2:nn-1)',uy(1:2:nn-1,1:2:nn-1)','w')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hGLL dhGLLdx] = LagrangeVal(xi,N,1);
eGLL = EdgeVal(dhGLLdx);

qeta = eGLL'*reshape(q(N*(N+1)+1:2*N*(N+1)),N+1,N)'*hGLL;
qxi  = hGLL'*reshape(q(1:N*(N+1)),N+1,N)*eGLL;

qx = (dxdxi.*qxi+dxdeta.*qeta)./J;
qy = (dydxi.*qxi+dydeta.*qeta)./J;

figure(10)
subplot(1,3,1)
pcolor(X,Y,qx)
shading interp
title('q_\xi')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

subplot(1,3,2)
pcolor(X,Y,qy)
shading interp
title('q_\eta')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

qq = sqrt(qx.^2+qy.^2);

subplot(1,3,3)
pcolor(X,Y,qq)
shading interp
title('q')
axis('square')
% colorbar
set(gca,'clim',[0 pi])
hold on
quiver(X(1:2:nn-1,1:2:nn-1),Y(1:2:nn-1,1:2:nn-1),qx(1:2:nn-1,1:2:nn-1)',qy(1:2:nn-1,1:2:nn-1)','w')