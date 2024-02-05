function PlotGlobalNumberingZeroFormDual(n,p,x,y,gnFlag, cornerFlag)
%PlotGlobalNumberingOneFormDual plots the Global numbering of a mesh of
%n(1) x n(2) elements of order p, given a global numbering gn.
%
%   PlotGlobalNumberingZeroFormDual(n,p,x,y,gnFlag)
%
%   n   :: is the number of elements. If n is a number then n x n elements
%          are assumed. If n is a vector with two entries then n(1) x n(2)
%          elements are assumed.
%   p   :: is the order of the mesh.
%   x   :: contains the x boundaries of the mesh
%   y   :: contains the y boundaries of the mesh
%   gnFlag  :: is a flag that states if global numbering is to be plotted or
%          not
%   cornerFlag  :: is a flag that states if corner nodes are to be plotted
%                  or not
%
%   See also: GLOBALNUMBERINGZEROFORMDUAL

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
    if gnFlag
        gn = GlobalNumberingZeroFormDual(n,p);
    end
    
    % compute the number of elements
    nElements = n(1)*n(2);
    
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
    
    % open a new figure
    figure()
    
    % set hold to on
    hold on
    
    if cornerFlag==false
        % loop over the elements and plot all the nodes and edges
        for k=1:nElements
            % scale the nodes
            x = 0.5*(xi+1)*deltax + lowerLeftNodes(k,1);
            y = 0.5*(eta+1)*deltay + lowerLeftNodes(k,2);

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

            if gnFlag
                % plot the interior node numbers
                text(xinterior(:), yinterior(:), num2cell(gn(k,1:(p*p))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center');
                % plot the node edges numbers
                text(x(2:(p+1),1),y(2:(p+1),1), num2cell(gn(k,(p*p+1):2:(p*p+2*p))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % left nodes

                text(x(2:(p+1),end),y(2:(p+1),end), num2cell(gn(k,(p*p+2):2:(p*p+2*p))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % right nodes
                
                text(x(1,2:(p+1)),y(1,2:(p+1)), num2cell(gn(k,(p*p+2*p+1):2:(end-4))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % bottom nodes
                
                text(x(end,2:(p+1)),y(end,2:(p+1)), num2cell(gn(k,(p*p+2*p+2):2:(end-4))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % top nodes
            end
        end
    else
        % loop over the elements and plot all the nodes and edges
        for k=1:nElements
            % scale the nodes
            x = 0.5*(xi+1)*deltax + lowerLeftNodes(k,1);
            y = 0.5*(eta+1)*deltay + lowerLeftNodes(k,2);

            % generate the element grid
            [x y] = meshgrid(x,y);

            % the interior nodes
            xinterior = x(2:(end-1),2:(end-1));
            yinterior = y(2:(end-1),2:(end-1));

            % plot all the nodes
            plot(x(:),y(:),'k.');

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

            if gnFlag
                % plot all the node numbers
                text(xinterior(:), yinterior(:), num2cell(gn(k,1:(p*p))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center');
                % plot the node edges numbers
                text(x(2:(p+1),1),y(2:(p+1),1), num2cell(gn(k,(p*p+1):2:(p*p+2*p))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % left nodes

                text(x(2:(p+1),end),y(2:(p+1),end), num2cell(gn(k,(p*p+2):2:(p*p+2*p))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % right nodes

                text(x(1,2:(p+1)),y(1,2:(p+1)), num2cell(gn(k,(p*p+2*p+1):2:(end-4))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % bottom nodes

                text(x(end,2:(p+1)),y(end,2:(p+1)), num2cell(gn(k,(p*p+2*p+2):2:(end-4))),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % top nodes
                
                text(reshape(x([1 end],[1 end]),[4 1]),reshape(y([1 end],[1 end]),[4 1]), num2cell(gn(k,(end-3):end)),'FontWeight','Bold','FontSize',12,...
                                                           'VerticalAlignment','middle',...
                                                           'HorizontalAlignment','center'); % top nodes
                
            end
        end
    end
end