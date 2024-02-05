% close all
if plot_fig || avimatlab

XYaxis = [ min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y)) ];

if ~exist('numElements','var'); numElements=1; end

for i=1:numElements

    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
%     surf(Xp,Yp,reshape(Up(:,i),nn,nn))
    pcolor(Xp,Yp,reshape(Up(:,i),nn,nn))
    hold on

end
for i=1:numElements
    Xfig = reshape(Mesh.X(:,i),N+1,N+1);
    Yfig = reshape(Mesh.Y(:,i),N+1,N+1);
    h1=mesh(Xfig,Yfig,0.6*ones(N+1)); %view([0 0 1]); axis equal
    set(h1,'FaceColor','None')
%     set(h1,'EdgeColor','None')
    hold on; grid off
end
shading interp
colorbar
axis equal
% axis([-1 1 -1 1 -0.2 1.2])
axis([-1 1 -1 1])
axis off
box off
% set(gca,'Ztick',[0 1])
set(gca,'clim',[0 1])
pause(0.05)
hold off
end

if avimatlab
    avimovie(filename,fig,t==0,t>=T);
end
% keyboard
if Tecplot
    data = [ dataXY Up ];
    if j<10
        name = strcat(filename,'_0',num2str(j));
    else
        name = strcat(filename,'_',num2str(j));
    end
    MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
end