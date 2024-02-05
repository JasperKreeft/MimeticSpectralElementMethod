function PlotMesh(phi,p,meshFigure,varargin)

    nElements = size(phi,1);
    
    xi = LobattoQuad(p+5);
    eta = xi;
    
    [Xi,Eta] = meshgrid(xi,eta);
    
    for element = 1:nElements
        
        [x,y] = phi{element}(Xi,Eta);
        z = ones(size(x));
        
        figure(meshFigure)
        if (size(varargin,2))
            plot3(x,y,varargin{1}*ones(size(x)),'-k','LineWidth',2)
            hold on
            plot3(x',y',varargin{1}*ones(size(x')),'-k','LineWidth',2)
        else
            plot(x,y,'-k','LineWidth',2)
            hold on
            plot(x',y','-k','LineWidth',2)
        end
        
    end
%     axis equal

end