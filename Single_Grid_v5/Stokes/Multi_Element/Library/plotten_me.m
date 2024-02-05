global nn

if isempty(nn)
    nn = sqrt(size(Meshp.X,1));
end

XYlim = [min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y))];
    
h1 = figure('visible','off');
% h2 = figure('visible','off');
% h3 = figure('visible','off');
% h4 = figure('visible','off');

for i=1:numElements
% 
    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
%     
    up = reshape(uu(:,i),nn,nn);
    vp = reshape(vv(:,i),nn,nn);
    Vp = reshape(velo(:,i),nn,nn);
    wp = reshape(ww(:,i),nn,nn);
    Pp = reshape(pp(:,i),nn,nn);
% %     divup = reshape(divu(:,i),nn,nn);
%     curlw1p = reshape(curlw1(:,i),nn,nn);
%     curlw2p = reshape(curlw2(:,i),nn,nn);
% 
    figure(h1)
    set(h1,'visible','off')
    subplot(3,2,1)
    pcolor(Xp,Yp,up)
    hold on
    shading interp
    % colorbar
    title('u')
    axis equal
    axis(XYlim)
    subplot(3,2,2)
    pcolor(Xp,Yp,vp)
    hold on
    shading interp
    % colorbar
    title('v')
    axis equal
    axis(XYlim)
    subplot(3,2,3)
    pcolor(Xp,Yp,Vp)
    hold on
    shading interp
    % colorbar
    title('abs(velo)')
    axis equal
    axis(XYlim)
    subplot(3,2,4)
    pcolor(Xp,Yp,wp)
    hold on
    shading interp
    % set(gca,'clim',[-5 5])
    % colorbar
    title('w')
    axis equal
    axis(XYlim)
    subplot(3,2,5)
    pcolor(Xp,Yp,Pp)
    hold on
    % set(gca,'clim',[-5 5])
    shading interp
    % colorbar
    title('p')
    axis equal
    axis(XYlim)
    subplot(3,2,6)
    mesh(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),zeros(N+1),'EdgeColor','black')
    hold on
    view([0 0 1])
    axis equal
    axis(XYlim)

%     figure(h2)
%     set(h2,'visible','off')
% %     pcolor(Xp,Yp,Vp)
% surf(Xp,Yp,Vp)
%     hold on
% %     quiver(Xp(1:4:nn-1,1:4:nn-1),Yp(1:4:nn-1,1:4:nn-1),up(1:4:nn-1,1:4:nn-1),vp(1:4:nn-1,1:4:nn-1),'w')
%     shading interp
% %     axis equal
% %     axis(XYlim)
% %     axis off
% %     set(gca,'clim',[0 1])
%     colorbar
%     title('Lid-Driven Cavity Stokes')

%     figure(h3)
%     set(h3,'visible','off')
%     surf(Xp,Yp,divup)
%     hold on
%     shading interp
%     colorbar
%     title('divergence u')
    
%     figure(h4)
%     set(h4,'visible','off')
%     subplot(1,2,1)
%     quiver(Xp,Yp,curlw1p,curlw2p)
% %     surf(Xp,Yp,curlw2p)
%     hold on
%     subplot(1,2,2)
%     curlw1exp = reshape(curlw1_ex(:,i),nn,nn);
%     curlw2exp = reshape(curlw2_ex(:,i),nn,nn);
%     quiver(Xp,Yp,curlw1exp-curlw1p,curlw2exp-curlw2p)
%     title('curl w')
%     hold on

end

set(h1,'visible','on')
% set(h2,'visible','on')
% set(h3,'visible','on')
% set(h4,'visible','on')
