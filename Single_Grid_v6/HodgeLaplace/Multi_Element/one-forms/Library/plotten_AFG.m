global nn

if isempty(nn)
    nn = sqrt(size(Meshp.X,1));
end

XYlim = [min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y))];
    
h1 = figure('visible','off');

for i=1:numElements

    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
    
    up = reshape(uu(:,i),nn,nn);
    vp = reshape(vv(:,i),nn,nn);
    divup = reshape(div_uv(:,i),nn,nn);
    sp = reshape(ssigma(:,i),nn,nn);
    ngsxp = reshape(ng_sx(:,i),nn,nn);%ng_x_ex
    ngsyp = reshape(ng_sy(:,i),nn,nn);%ng_y_ex

%     divup = reshape(divu(:,i),nn,nn);
%     curlw1p = reshape(curlw1(:,i),nn,nn);
%     curlw2p = reshape(curlw2(:,i),nn,nn);

    figure(h1)
    set(h1,'visible','off')
    subplot(2,3,1)
    pcolor(Xp,Yp,up)
    hold on
    shading interp
%     colorbar
    title('u')
    axis equal
    axis(XYlim)
    subplot(2,3,2)
    pcolor(Xp,Yp,vp)
    hold on
    shading interp
%     colorbar
    title('v')
    axis equal
    axis(XYlim)
    subplot(2,3,3)
    pcolor(Xp,Yp,divup)
    hold on
    shading interp
    % colorbar
    title('div u')
    axis equal
    axis(XYlim)
    subplot(2,3,4)
    pcolor(Xp,Yp,sp)
    hold on
    shading interp
    % set(gca,'clim',[-5 5])
    % colorbar
    title('\sigma')
    axis equal
    axis(XYlim)
    subplot(2,3,5)
    pcolor(Xp,Yp,ngsxp)
    hold on
    % set(gca,'clim',[-5 5])
    shading interp
    % colorbar
    title('(d\sigma)_x')
    axis equal
    axis(XYlim)
    subplot(2,3,6)
    pcolor(Xp,Yp,ngsyp)
    hold on
    % set(gca,'clim',[-5 5])
    shading interp
    % colorbar
    title('(d\sigma)_y')
    axis equal
    axis(XYlim)    
%     mesh(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),zeros(N+1),'EdgeColor','black')
%     hold on
%     view([0 0 1])
%     axis equal
%     axis(XYlim)

end
    
set(h1,'visible','on')
