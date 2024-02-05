global c m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% postproces coordinates and integration weights 
nn = 100;
[xi,wg] = Gnodes(nn); eta = xi;
Xi = xi'*ones(1,nn);  Eta = Xi';
Wgg = wg'*wg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hGL eGL] = MimeticpolyVal(xi,N,1);
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

phi = reshape(phi_in,N,N);

PHI = reshape(phi,[],1);

figure(9)
subplot(2,3,1)
surf(XGG,YGG,phi)
view([0 0 1])
axis('square')
title('\phi')

Qxi  = reshape(q(1:N*(N+1)),N+1,N);
Qeta = reshape(q(N*(N+1)+1:2*N*(N+1)),N+1,N)';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pphi = (eGL'*phi*eGL).*J;

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

qxi  =  hGL'*Qxi*eGL;
qeta =  eGL'*Qeta*hGL;
qx = (qxi.*dxdxi+qeta.*dxdeta)./J;
qy = (qxi.*dydxi+qeta.*dydeta)./J;

figure(9)
subplot(2,3,4)
pcolor(X,Y,qx)
shading interp
title('u')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

figure(9)
subplot(2,3,5)
pcolor(X,Y,qy)
shading interp
title('v')
axis('square')
% colorbar
set(gca,'clim',[-pi pi])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

velo = sqrt(qx.^2+qy.^2);

figure(9)
subplot(2,3,6)
pcolor(X,Y,velo)
shading interp
axis('square')
% colorbar
set(gca,'clim',[0 pi])
title('absolute velocity')
hold on
quiver(X(1:2:nn-1,1:2:nn-1),Y(1:2:nn-1,1:2:nn-1),qx(1:2:nn-1,1:2:nn-1)',qy(1:2:nn-1,1:2:nn-1)','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%