function PlotGlobalNumberingOneFormTwoD(n,p,x,y,gnFlag,active)
%
%   PlotGlobalNumberingOneFormTwoD plots the Global numbering of a mesh of
%   n(1) x n(2) elements of order p, given a global numbering gn.
%
%   PlotGlobalNumberingOneFormTwoD(n,p,x,y,gnFlag)
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
    if any(gnFlag)
        gn = GlobalNumberingOneFormTwoD(n,p);
    end
    
    % compute the number of elements
    nElements = n(1)*n(2);
    
    % compute the number of edges in each element
    nElementEdgesX = p(1)*(p(2)+1);
    % nElementEdgesY = p(2)*(p(1)+1);
    
    % compute domain size
    Lx = x(2)-x(1);
    Ly = y(2)-y(1);
    
    % element lengths
    deltax = Lx/n(1); 
    deltay = Ly/n(2); 
    
    % allocate memory space for lower left node of elements
    lowerLeftNodes = zeros(n(1)*n(2),2);
    
    % compute the x coordinates of the lower left nodes
    lowerLeftNodes(:,1) = repmat(x(1)+(0:(n(1)-1))*deltax,[1 n(2)])';
    
    % compute the y coordinates of the lower left nodes
    lowerLeftNodes(:,2) = rectpulse(y(1)+(0:(n(2)-1))*deltay,n(1))';
    
    % generate element nodes
    xi = GLLnodes(p(1));
    eta = GLLnodes(p(2));
    
    % generate edge coordinates (xi and eta)
    xiEdgedXi = GLLnodes(p(1));
    xiEdgedXi = 0.5*(xiEdgedXi(2:end)-xiEdgedXi(1:p(1)))+xiEdgedXi(1:p(1));
    etaEdgedXi = GLLnodes(p(2));
    
    xiEdgedEta = GLLnodes(p(1));
    etaEdgedEta = GLLnodes(p(2));
    etaEdgedEta = 0.5*(etaEdgedEta(2:end)-etaEdgedEta(1:p(2)))+etaEdgedEta(1:p(2));
    
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
        
        % scale the edge text location
        xEdgedXi = 0.5*(xiEdgedXi+1)*deltax+lowerLeftNodes(k,1);
        yEdgedXi = 0.5*(etaEdgedXi+1)*deltay+lowerLeftNodes(k,2);        
        xEdgedEta = 0.5*(xiEdgedEta+1)*deltax+lowerLeftNodes(k,1);
        yEdgedEta = 0.5*(etaEdgedEta+1)*deltay+lowerLeftNodes(k,2);
        
        % generate the element grid
        [x y] = meshgrid(x,y);
        
        % generate the edge text locations grid
        [xEdgedXi yEdgedXi] = meshgrid(xEdgedXi,yEdgedXi);
        [xEdgedEta yEdgedEta] = meshgrid(xEdgedEta,yEdgedEta);
        
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
        
        if gnFlag(1)
            % plot the edge numbers dXi
            text(reshape(xEdgedXi',numel(xEdgedXi),1),reshape(yEdgedXi',numel(yEdgedXi),1),num2cell(gn(k,1:nElementEdgesX)),...
                 'FontWeight','Bold','FontSize',12,'VerticalAlignment','bottom','HorizontalAlignment','center','color','red');
        end
        
        if gnFlag(2)
            % plot the edge numbers dEta
            text(xEdgedEta(:),yEdgedEta(:),num2cell(gn(k,(nElementEdgesX+1):end)),...
                 'FontWeight','Bold','FontSize',12,'VerticalAlignment','middle','HorizontalAlignment','left','color','green');
        end
        
    end
    
    % extra plot options
    axis equal; axis tight; 
    % set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);

end