function PlotGlobalNumberingTwoFormTwoD(n,p,x,y,gnFlag,active)
%
%   PlotGlobalNumberingTwoFormTwoD plots the Global numbering of a mesh of
%   n(1) x n(2) elements of order p, given a global numbering gn.
%
%   PlotGlobalNumberingTwoFormTwoD(n,p,x,y,gnFlag)
%
%   input:
%       n       :: the number of elements; if n is a number it has n x n
%                  elements, if it is a vector, then it has n(1) x n(2) elements
%       p       :: the order of the 0-form approximation; if p is a number then 
%                  px=py=p, if it is a vector, then px!=py
%       x       :: contains the x boundaries of the mesh
%       y       :: contains the y boundaries of the mesh
%       gnFlag  :: is a flag that states if global numbering is to be plotted or not
%       active  :: flag for using active figure or not
%
%   output:
%       ...     :: plot
%
%   Copyright 2010 Artur Palha
%   $Revision: 1.2 $  $Date: 02/03/2010 $
%   Revised by: Peter Kuystermans
%   $Revision: 2.0 $  $Date: 13/07/2011 $ 
%   Revised by: Peter Kuystermans
%   $Revision: 3.0 $  $Date: 17/10/2011 $ 
    
%-------------------------------------------------------------------------%
% input checks                                                            %
%-------------------------------------------------------------------------%       
 
    if length(n)>2
        fprintf(':: n :: has too many values \n');
        return
    elseif length(n)==1
        n = [n n];
    end   
    if prod(n)==0
        fprintf(':: n :: has at least one zero \n');
        return
    end
    
    if length(p)>2
        fprintf(':: p :: has too many values \n');
        return
    elseif length(p)==1
        p = [p p];
    end    
    if prod(p)==0
        fprintf(':: p :: has at least one zero \n');
        return
    end      
    
%-------------------------------------------------------------------------%
% calculation of coordinates                                              %
%-------------------------------------------------------------------------%
    
    % generate the global numbering if it is to be plotted
    if gnFlag
        gn = GlobalNumberingTwoFormTwoD(n,p);
    end
    
    % compute the number of elements
    nElements = n(1)*n(2);
    
    % compute the number of surfaces in each element
    nElementSurfaces = p(1)*p(2);
    
    % compute domain size
    Lx = x(2)-x(1);
    Ly = y(2)-y(1);
    
    % element lengths
    deltax = Lx/n(1); 
    deltay = Ly/n(2); 
    
    % allocate memory space for lower left node of elements
    lowerLeftNodes = zeros(n(1)*n(2),2);
    
    % compute the x coordinates of the lower left nodes
    lowerLeftNodes(:,1) = repmat(x(1)+(0:(n(1)-1))*deltax, [1 n(2)])';
    
    % compute the y coordinates of the lower left nodes
    lowerLeftNodes(:,2) = rectpulse(y(1)+(0:(n(2)-1))*deltay, n(1))';
    
    % generate element nodes
    xi = GLLnodes(p(1));
    eta = GLLnodes(p(2));
    
    % calculation of surface centers for standard element
    xStandard = GLLnodes(p(1));
    xStandard = 0.5*(xStandard(2:end)-xStandard(1:p(1)))+xStandard(1:p(1));
    yStandard = GLLnodes(p(2));
    yStandard = 0.5*(yStandard(2:end)-yStandard(1:p(2)))+yStandard(1:p(2));
    
%-------------------------------------------------------------------------%
% the plotting                                                            %
%-------------------------------------------------------------------------%  
    
    % open a new figure
    if active == 1
        current_figure = gcf;
        figure(current_figure)
    else
        figure()
    end
    
    % set hold to on
    hold on
    
    % loop over the elements and plot all the nodes and edges
    for k=1:nElements
        
        % scale the nodes
        x = 0.5*(xi+1)*deltax+lowerLeftNodes(k,1);
        y = 0.5*(eta+1)*deltay+lowerLeftNodes(k,2);
        
        % generate the element grid
        [x y] = meshgrid(x,y);
        
        % scale the surface text location
        xSurf = 0.5*(xStandard+1)*deltax+lowerLeftNodes(k,1);
        ySurf = 0.5*(yStandard+1)*deltay+lowerLeftNodes(k,2);
        
        % generate the surface text grid
        [xSurf ySurf] = meshgrid(xSurf,ySurf);
        
        % plot the nodes
        plot(x(:),y(:),'k.');
        
        % plot the interior edges
        for n=2:p(2)
            plot(x(n,:),y(n,:),'k-');
        end
        for n=2:p(1)
            plot(x(:,n),y(:,n),'k-');
        end

        % plot the boundary edges
        for n=[1 p(2)+1]
            plot(x(n,:),y(n,:),'b-','LineWidth',2);
        end
        for n=[1 p(1)+1]
        	plot(x(:,n),y(:,n),'b-','LineWidth',2);
        end
        
        if gnFlag
            text(reshape(xSurf',numel(xSurf),1),reshape(ySurf',numel(ySurf),1),num2cell(gn(k,1:nElementSurfaces)),...
                 'FontWeight','Bold','FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','center','color','red');
        end
        
    end
    
    % extra plot options
    axis equal; axis tight; 
    % set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);

end