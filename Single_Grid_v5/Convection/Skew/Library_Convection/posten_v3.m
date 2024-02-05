global nn

if plot_fig || avimatlab

XYaxis = [ min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y)) ];

for i=1:numElements

    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
%     surf(Xp,Yp,reshape(Up(:,i),nn,nn))
    pcolor(Xp,Yp,reshape(Up(:,i),nn,nn))
    hold on

end
shading interp
colorbar
axis equal
% axis([-1 1 -1 1 -0.2 1.2])
axis(XYaxis)
% set(gca,'Ztick',[0 1])
set(gca,'clim',[0 1])
pause(0.05)
hold off
end

if avimatlab
    avimovie(filename,fig,t==0,t>=T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Tecplot

    data = zeros(nn^2*numElements,3);
    corners = zeros((nn-1)^2*numElements,4);
    for nel=1:numElements

        for i=1:nn-1
            for k=1:nn-1

                iknel = i+(k-1)*(nn-1)+(nel-1)*(nn-1)^2;

                ll = i+(k-1)*nn; % lowerleft corner

                corners(iknel,:) = [ ll ll+1 ll+nn+1 ll+nn ] + (nel-1)*nn^2;

            end
        end


        ind = (1:nn^2)+(nel-1)*nn^2;
        data(ind,:) = [ Meshp.X(:,nel) Meshp.Y(:,nel) Up(:,nel) ];

        if j<10
            name = strcat(filename,'_0',num2str(j));
        else
            name = strcat(filename,'_',num2str(j));
        end

    end

    MatlabToTecplot('FE',name,name,'"X" "Y" "U"',[numElements*nn^2 numElements*(nn-1)^2],data,2,corners);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*IJCL>