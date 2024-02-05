
if plot_fig || avimatlab
% surf(Xp,Yp,reshape(Up,nn,nn))
pcolor(Xp,Yp,reshape(Up,nn,nn))
shading interp
colorbar
axis equal
% axis([-1 1 -1 1 -0.2 1.2])
axis([-1 1 -1 1])
% set(gca,'Ztick',[0 1])
set(gca,'clim',[0 1])
pause(0.05)
end

if avimatlab
    avimovie(filename,fig,t==dt,t>=T);
end

if Tecplot
    data = [ dataXY reshape(Up,[],1) ];
    if j<10
        name = strcat(filename,'_0',num2str(j));
    else
        name = strcat(filename,'_',num2str(j));
    end
    MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
end