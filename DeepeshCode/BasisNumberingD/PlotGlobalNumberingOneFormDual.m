function PlotGlobalNumberingOneFormDual(n,p,x,y,gnTypeFlag,gnPlotFlag)
%PlotGlobalNumberingOneFormDual plots the Global numbering of a mesh of
%n(1) x n(2) elements of order p, given a global numbering gn.
%
%   PlotGlobalNumberingOneFormDual(n,p,x,y,gnPlotFlag)
%
%   n   :: is the number of elements. If n is a number then n x n elements
%          are assumed. If n is a vector with two entries then n(1) x n(2)
%          elements are assumed.
%   p   :: is the order of the mesh.
%   x   :: contains the x boundaries of the mesh
%   y   :: contains the y boundaries of the mesh
%   gnTypeFlag  :: is a flag that states if global numbering is to be done
%                  such that adjacent eges of two elements should share the
%                  same number (true) or not (false)
%   gnPlotFlag  :: is a flag that states if global numbering is to be plotted or
%          not
%
%   See also: GLOBALNUMBERINGONEFORMDUAL, GLOBALNUMBERINGZEROFORMDUAL

%   Copyright 2010 Artur Palha
%   $Revision: 1.2 $  $Date: 2010/03/08 $
    
    if length(n)>2
        % check if n has more than two values
        disp(sprintf(':: n :: has too many values'));
        return
    elseif length(n)==1
            n = [n n];
    end
    
    if prod(n)==0
        % check if n has a zero value
        disp(sprintf(':: n :: has at least one zero'));
        return
    end
    
    % generate the global numbering if it is to be plotted
    if any(gnPlotFlag)
        if gnTypeFlag
            gn = GlobalNumberingOneFormDualAlt(n,p);
        else
            gn = GlobalNumberingOneFormDual(n,p);
        end
    end
    
    % compute the number of elements
    nElements = n(1)*n(2);
    
    % compute the number of edges in each element
    nElementEdges = p*(p+1);
    
    % compute the location of the corner nodes of each element
    Lx = x(2)-x(1);
    Ly = y(2)-y(1);
    
    deltax = Lx/n(1); % x length of elements
    deltay = Ly/n(2); % y length of elements
    
    % allocate memory space for lower left node of elements
    lowerLeftNodes = zeros(n(1)*n(2),2);
    
    % compute the x coordinates of the lower left nodes
    lowerLeftNodes(:,1) = rectpulse(x(1)+(0:(n(1)-1))*deltax, n(2))';
    
    % compute the y coordinates of the lower left nodes
    lowerLeftNodes(:,2) = repmat(y(1)+(0:(n(2)-1))*deltay, [1 n(1)])';
    
    % generate element nodes in canonical intervals
    xi = EGaussQuad(p+1);
    eta = xi;
    
    % generate edge coordinates in canonical intervals
    xiEdgedXi = EGaussQuad(p+1);
    xiEdgedXi = 0.5*(xiEdgedXi(2:end)-xiEdgedXi(1:(p+1))) + xiEdgedXi(1:(p+1));
    etaEdgedXi = xi;
    
    xiEdgedEta = etaEdgedXi;
    etaEdgedEta = xiEdgedXi;
    
    % open a new figure
    figure()
    
    % set hold to on
    hold on
        
    % loop over the elements and plot all the nodes and edges
    for k=1:nElements
        % scale the nodes
        x = 0.5*(xi+1)*deltax + lowerLeftNodes(k,1);
        y = 0.5*(eta+1)*deltay + lowerLeftNodes(k,2);
        
        % scale the edge text location
        xEdgedXi = 0.5*(xiEdgedXi+1)*deltax + lowerLeftNodes(k,1);
        yEdgedXi = 0.5*(etaEdgedXi+1)*deltay + lowerLeftNodes(k,2);
        
        xEdgedEta = 0.5*(xiEdgedEta+1)*deltax + lowerLeftNodes(k,1);
        yEdgedEta = 0.5*(etaEdgedEta+1)*deltay + lowerLeftNodes(k,2);
        
        % generate the element grid
        [x y] = meshgrid(x,y);
        
        % the interior nodes
        xinterior = x(2:(end-1),2:(end-1));
        yinterior = y(2:(end-1),2:(end-1));
        
        % plot the interior nodes
        plot(xinterior(:),yinterior(:),'k.');
        
        % plot the boundary nodes
        plot(x(1,2:(p+1)),y(1,2:(p+1)),'k.'); % bottom
        plot(x(end,2:(p+1)),y(end,2:(p+1)),'k.'); % top
        plot(x(2:(p+1),1),y(2:(p+1),1),'k.'); % left
        plot(x(2:(p+1),end),y(2:(p+1),end),'k.'); % right
        
        % plot the interior edges
        for n=2:(p+1)
            plot(x(n,:),y(n,:),'k-');
            plot(x(:,n),y(:,n),'k-');
        end
        
        % plot the boundary edges
        for n=[1 p+2]
            plot(x(n,:),y(n,:),'k-','LineWidth',1.5);
            plot(x(:,n),y(:,n),'k-','LineWidth',1.5);
        end
        
        if gnPlotFlag(1) % xi edges
            %generate the edge text locations grid
            [xEdgedXi yEdgedXi] = meshgrid(xEdgedXi,yEdgedXi(2:(p+1)));
            
            % plot the edge numbers dXi
            text(xEdgedXi(:), yEdgedXi(:), num2cell(gn(k,(nElementEdges+1):end)),...
                                          'FontWeight','Bold','FontSize',12,...
                                          'VerticalAlignment','bottom',...
                                          'HorizontalAlignment','center');
        end
        
        if gnPlotFlag(2) % eta edges
            [xEdgedEta yEdgedEta] = meshgrid(xEdgedEta(2:(p+1)),yEdgedEta);
            % plot the edge numbers dEta
            text(xEdgedEta(:), yEdgedEta(:), num2cell(gn(k,1:nElementEdges)),...
                                            'FontWeight','Bold','FontSize',12,...
                                            'VerticalAlignment','middle',...
                                            'HorizontalAlignment','left');
        end
    end

end