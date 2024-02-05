XYaxis = [ min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y)) ];

figure
% keyboard
nn = sqrt(length(Meshp.X(:,1)));

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
