global nn

XYaxis = [ min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y)) ];

figure

for i=1:numElements

    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
    
    if exist('phi_ex','var')
    subplot(2,3,1)
    phip = reshape(phi_ex(:,i),nn,nn);
    surf(Xp,Yp,phip)
    hold on
    end

    subplot(2,3,2)
    phip = reshape(phi(:,i),nn,nn);
    surf(Xp,Yp,phip)
    hold on

    subplot(2,3,3)
    if exist('phi_interp','var')
    phip = reshape(phi_interp(:,i),nn,nn);
    surf(Xp,Yp,phip)
    hold on
    end

    if exist('phi_ex','var')
    subplot(2,3,4)
    phip = reshape(phi_ex(:,i)-phi(:,i),nn,nn);
    surf(Xp,Yp,phip)
    hold on
    end

    if exist('phi_interp','var')
    subplot(2,3,5)
    phip = reshape(phi_interp(:,i)-phi(:,i),nn,nn);
    surf(Xp,Yp,phip)
    hold on
    end
    
    if exist('phi_ex','var') && exist('phi_interp','var')
    subplot(2,3,6)
    phip = reshape(phi_ex(:,i)-phi_interp(:,i),nn,nn);
    surf(Xp,Yp,phip)
    hold on
    end

end

subplot(2,3,1)
view([0 0 1])
title('\phi_{ex}')
% axis([-1 1 -1 1])
% axis('square')
shading interp
% colorbar
% set(gca,'clim',[-1 1])
% set(gca,'clim',[0 1])

subplot(2,3,2)
view([0 0 1])
title('\phi')
% axis([-1 1 -1 1])
% axis('square')
shading interp
% colorbar
% set(gca,'clim',[-1 1])

subplot(2,3,3)
view([0 0 1])
title('\phi_{interp}')
% axis([-1 1 -1 1])
% axis('square')
shading interp
% colorbar
% set(gca,'clim',[-1 1])
% set(gca,'clim',[0 1])

subplot(2,3,4)
view([0 0 1])
title('\phi_{ex}-\phi')
% axis([-1 1 -1 1])
% axis('square')
shading interp
colorbar

subplot(2,3,5)
view([0 0 1])
title('\phi_{interp}-\phi')
% axis([-1 1 -1 1])
% axis('square')
shading interp
colorbar

subplot(2,3,6)
view([0 0 1])
title('\phi_{ex}-\phi_{interp}')
% axis([-1 1 -1 1])
% axis('square')
shading interp
colorbar

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

for i=1:numElements

    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);

    subplot(2,3,1)
    Qp = reshape(qx(:,i),nn,nn);
    surf(Xp,Yp,Qp)
    hold on

    subplot(2,3,2)
    Qp = reshape(qy(:,i),nn,nn);
    surf(Xp,Yp,Qp)
    hold on

    subplot(2,3,3)
    Qp = reshape(qMag(:,i),nn,nn);
    surf(Xp,Yp,Qp)
    hold on

    if exist('qx_ex','var')
    subplot(2,3,4)
    Qp = reshape(qx_ex(:,i),nn,nn);
    surf(Xp,Yp,Qp)
    hold on

    subplot(2,3,5)
    Qp = reshape(qy_ex(:,i),nn,nn);
    surf(Xp,Yp,Qp)
    hold on
    end

    subplot(2,3,6)
    Qp = reshape(qMag(:,i),nn,nn);
    contourf(Xp,Yp,Qp)
    axis equal
    axis(XYaxis)
    hold on
end

subplot(2,3,1)
view([0 0 1])
title('q_x')
shading interp
colorbar

subplot(2,3,2)
view([0 0 1])
title('q_y')
shading interp
colorbar

subplot(2,3,3)
view([0 0 1])
title('Mag q')
shading interp
colorbar

subplot(2,3,4)
view([0 0 1])
title('q_x exact')
shading interp
colorbar

subplot(2,3,5)
view([0 0 1])
title('q_y exact')
shading interp
colorbar

subplot(2,3,6)
view([0 0 1])
title('empty')
axis equal
axis(XYaxis)
% shading interp
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(FunctionType,'nozzle') && strcmp(Domain,'LavalNozzle')
    figure
    axis equal
    axis([-1 1 0 1])
    hold on
    for i=1:numElements
    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
    Qpx = reshape(qx(:,i),nn,nn);
    Qpy = reshape(qy(:,i),nn,nn);
    Qp = reshape(qMag(:,i),nn,nn);
    pcolor(Xp,Yp,Qp)
    title('absolute flux')
    quiver(Xp(1:2:end,1:2:end),Yp(1:2:end,1:2:end),Qpx(1:2:end,1:2:end),Qpy(1:2:end,1:2:end),'w')
    colorbar
    end
    shading interp
end