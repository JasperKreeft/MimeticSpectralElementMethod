global nn

if plot_fig || avimatlab

if ~exist('numElements','var')
    numElements = 1;
end

XYaxis = [ min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y)) ];

Xp = reshape(Meshp.X(:,1),nn,nn);
Yp = reshape(Meshp.Y(:,1),nn,nn);
% surf(Xp,Yp,reshape(Wp(:,1),nn,nn))
pcolor(Xp,Yp,reshape(Wp(:,1),nn,nn))
hold on
for i=2:numElements

    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
%     surf(Xp,Yp,reshape(Wp(:,i),nn,nn))
    pcolor(Xp,Yp,reshape(Wp(:,i),nn,nn))

end
shading interp
colorbar
axis equal
axis(XYaxis)
% axis([-1 1 -1 1 -0.1 1.])
% set(gca,'Ztick',[0 1])
set(gca,'clim',Clim)
pause(0.05)
hold off
end


if avimatlab
    set(gcf,'visible','off')
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
        data(ind,:) = [ Meshp.X(:,nel) Meshp.Y(:,nel) Wp(:,nel) ];

        if j<10
            name = strcat(filename,'_0',num2str(j));
        else
            name = strcat(filename,'_',num2str(j));
        end

    end

    MatlabToTecplot('FE',name,name,'"X" "Y" "W"',[numElements*nn^2 numElements*(nn-1)^2],data,2,corners);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%