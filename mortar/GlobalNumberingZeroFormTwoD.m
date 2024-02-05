function gn = GlobalNumberingZeroFormTwoD(n,p)
%
%   GlobalNumberingZeroFormTwoD computes the global numbering of 0-form
%   degrees of freedom on a primal mesh.
%
%   gn = GlobalNumberingZeroFormTwoD(n,p)
%
%   input:
%       n   :: the number of elements; if n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements
%       p   :: the order of the 0-form approximation; if p is a number then 
%              px=py=p, if it is a vector, then px!=py
%
%   output:
%       gn  :: global numbering of the nodes of the mesh (gn is a N x (p+1)^2 
%              matrix where N is the total number of elements)
%
%   Local (per element) numbering:
%
%       1 - surfaces are numbered from left to right and bottom to top
%       2 - horizontal edges are numbered from left to right and bottom to top    
%       3 - vertical edges are numbered from bottom to top and left to right
%       4 - nodes are numbered from left to right and bottom to top
%
%   Global numbering:
%
%       note: within each element local numbering rules apply
%       elements:
%       1.1 - elements are numbered from left to right and bottom to top 
%       surfaces
%       2.1 - next the element surfaces are numbered
%       edges
%       3.1 - internal horizontal edges of the elements
%       3.2 - internal horizontal element boundary edges
%       3.3 - internal vertical edges of the elements
%       3.4 - internal vertical element boundary edges
%       3.5 - bottom/top domain boundaries
%       3.6 - left/right domain boundaries
%       nodes
%       4.1 - internal element nodes
%       4.2 - internal element corner nodes
%       4.3 - internal horizontal element boundary nodes
%       4.4 - internal vertical element boundary nodes
%       4.5 - bottom/top domain boundary nodes (corners excluded)
%       4.6 - left/right domain boundary nodes (corners excluded)
%       4.7 - (bottom/top) domain corner nodes
%       4.8 - bottom/top internal domain boundary nodes
%       4.9 - left/right internal domain boundary nodes
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
% parameter values                                                        %
%-------------------------------------------------------------------------%
    
    % compute the number of elements
    nElements = n(1)*n(2); % number of elements

    % compute the number of interior elements
    if prod(n)==1
        nInternalElements = 0;
    else
        nInternalElements = (n(1)-2)*(n(2)-2);
    end

    % compute the number of internal nodes
    nInternalNodes = (p(1)-1)*(p(2)-1)*nElements;

    % compute the number of internal nodes per element
    nElementInternalNodes = (p(1)-1)*(p(2)-1);

    % compute the number of nodes in each element
    nElementNodes = (p(1)+1)*(p(2)+1);

    % allocate memory space for global numbering matrix
    gn = zeros(nElements,nElementNodes,'uint32');

%-------------------------------------------------------------------------%
% number the internal element nodes                                       %
%-------------------------------------------------------------------------%

    % select the internal nodes
    internalNodesColumns = zeros(nElementInternalNodes,1,'uint32');
    count = 1;
    for i=2:p(2)
        for j=2:p(1)   
            internalNodesColumns((i-2)*(p(1)-1)+j-1) = (p(1)+1)*(i-1)+j;
            count = count+1;           
        end
    end
    
    % number the internal nodes
    if nInternalNodes==0
        nodesgn = 0;
    else
        nodesgn = 1:nInternalNodes;
        gn(:,internalNodesColumns) = reshape(nodesgn,[(p(1)-1)*(p(2)-1) nElements])';
    end

%-------------------------------------------------------------------------%
% number corner boundary nodes of the interior elements                   %
%-------------------------------------------------------------------------%   

    % compute the global numbering of the nodes
    if (n(1)-1)*(n(2)-1)>0
        nodesgn = uint32(nInternalNodes+1):uint32(nInternalNodes+(n(1)-1)*(n(2)-1));
    end

    % compute the internal elements
    internalElements = zeros(nInternalElements,1,'uint32');
    for i=2:(n(2)-1)
        for j=2:(n(1)-1)                    
            internalElements((i-2)*(n(1)-2)+j-1) = (i-1)*(n(1))+j;
        end 
    end 

    % compute the internal numbering of the corner nodes
    cornerNodesInternalNumber = [1 p(1)+1 p(2)*(p(1)+1)+1 (p(1)+1)*(p(2)+1)];

    % number the nodes on the inner elements
    if nInternalElements >= 1       
        for k=1:(n(2)-2)        
            gn(internalElements(((k-1)*(n(1)-2)+1):(k*(n(1)-2))),cornerNodesInternalNumber(1)) = nodesgn(((k-1)*(n(1)-1)+1):(k*(n(1)-1)-1));           
            gn(internalElements(((k-1)*(n(1)-2)+1):(k*(n(1)-2))),cornerNodesInternalNumber(2)) = nodesgn(((k-1)*(n(1)-1)+1+1):(k*(n(1)-1)-1+1));           
            gn(internalElements(((k-1)*(n(1)-2)+1):(k*(n(1)-2))),cornerNodesInternalNumber(3)) = nodesgn(((k)*(n(1)-1)+1):((k+1)*(n(1)-1)-1));
            gn(internalElements(((k-1)*(n(1)-2)+1):(k*(n(1)-2))),cornerNodesInternalNumber(4)) = nodesgn(((k)*(n(1)-1)+1+1):((k+1)*(n(1)-1)-1+1));  
        end     
    end

    % update boundary elements (excluding corner elements)
    if (n(2)>1)&&(n(1)>2)         
        % elements on the lower boundary
        gn(2:(n(1)-1),cornerNodesInternalNumber(3)) = nodesgn(1:(n(1)-2));
        gn(2:(n(1)-1),cornerNodesInternalNumber(4)) = nodesgn(2:(n(1)-1));
        % elements on the upper boundary
        gn(((n(2)-1)*n(1)+2):(n(1)*n(2)-1),cornerNodesInternalNumber(1)) = nodesgn(((n(1)-1)*(n(2)-2)+1):((n(1)-1)*(n(2)-1)-1));
        gn(((n(2)-1)*n(1)+2):(n(1)*n(2)-1),cornerNodesInternalNumber(2)) = nodesgn(((n(1)-1)*(n(2)-2)+2):((n(1)-1)*(n(2)-1)));
    end

    if (n(1)>1)&&(n(2)>2) 
        %elements on the left boundary
        gn((n(1)+1):n(1):((n(2)-1)*n(1)),cornerNodesInternalNumber(2)) = nodesgn((1:n(1)-1:((n(2)-2)*(n(1)-1))));
        gn((n(1)+1):n(1):((n(2)-1)*n(1)),cornerNodesInternalNumber(4)) = nodesgn((1:n(1)-1:((n(2)-2)*(n(1)-1)))+n(1)-1);
        % elements on the right boundary
        gn((2*n(1)):n(1):((n(2)-1)*n(1)),cornerNodesInternalNumber(1)) = nodesgn((1+n(1)-2):n(1)-1:((n(2)-2)*(n(1)-1)));
        gn((2*n(1)):n(1):((n(2)-1)*n(1)),cornerNodesInternalNumber(3)) = nodesgn(((1+n(1)-2):n(1)-1:((n(2)-2)*(n(1)-1)))+n(1)-1);
    end 
    
    % update corner elements
    if (n(1)-1)*(n(2)-1)>0        
        % left lower corner
        gn(1,cornerNodesInternalNumber(4)) = nodesgn(1); 
        if n(1)>1
            % right lower corner
            gn(n(1),cornerNodesInternalNumber(3)) = nodesgn(n(1)-1); 
        end
        if n(2)>1
            % left upper corner
            gn((n(2)-1)*n(1)+1,cornerNodesInternalNumber(2)) = nodesgn((n(2)-2)*(n(1)-1)+1);
        end
        % right upper corner
        gn(n(1)*n(2),cornerNodesInternalNumber(1)) = nodesgn(end);    
    end

%-------------------------------------------------------------------------%
% number the nodes of the internal edges                                  %
%-------------------------------------------------------------------------%

    leftEdgeNodesInternalNumber = (p(1)+2):p(1)+1:(p(2)-1)*(p(1)+1)+1; 
    rightEdgeNodesInternalNumber = leftEdgeNodesInternalNumber+p(1);
    bottomEdgeNodesInternalNumber = 2:p(1);
    topEdgeNodesInternalNumber = (p(2)*(p(1)+1)+2):((p(1)+1)*(p(2)+1)-1);
    
    if p(1)>1
        % horizontal internal boundary nodes
        if n(2)~=1
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p(1)-1)*n(1)*(n(2)-1));
            for k=1:(n(2)-1)
                gn(((k-1)*n(1)+1):k*n(1),topEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p(1)-1)*n(1)+1):(k*n(1)*(p(1)-1))),[p(1)-1 n(1)])';
                gn((k*n(1)+1):(k+1)*n(1),bottomEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p(1)-1)*n(1)+1):(k*n(1)*(p(1)-1))),[p(1)-1 n(1)])';           
            end
        end
    end
    if p(2)>1
        % vertical internal boundary nodes
        if n(1)~=1
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p(2)-1)*n(2)*(n(1)-1));
            for k=1:(n(1)-1)
                gn(k:n(1):(n(1)*n(2)),rightEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p(2)-1)*n(2)+1):(k*n(2)*(p(2)-1))),[p(2)-1 n(2)])';
                gn((k:n(1):(n(1)*n(2)))+1,leftEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p(2)-1)*n(2)+1):(k*n(2)*(p(2)-1))),[p(2)-1 n(2)])';
            end
        end
    end

%-------------------------------------------------------------------------%
% number the nodes at the domain boundary                                 %
%-------------------------------------------------------------------------%

    lowerLeft = 1;
    upperLeft = p(2)*(p(1)+1)+1;
    lowerRight = p(1)+1;
    upperRight = (p(1)+1)*(p(2)+1);

    % number the element corner nodes (domain corner nodes excluded)
    if n(1)>1
        % lower boundary
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)-1);
        gn(1:(n(1)-1),lowerRight) = nodesgn'; 
        gn(2:n(1),lowerLeft) = nodesgn'; 
        % upper boundary  
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)-1);
        gn(((n(2)-1)*n(1)+1):(n(1)*n(2)-1),upperRight) = nodesgn'; 
        gn(((n(2)-1)*n(1)+2):(n(1)*n(2)),upperLeft) = nodesgn'; 
    end    
    if n(2)>1
        % left boundary
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)-1);
        gn(1:n(1):((n(2)-1)*n(1)),upperLeft) = nodesgn'; 
        gn((n(1)+1):n(1):(n(1)*n(2)),lowerLeft) = nodesgn'; 
        % right boundary
        nodesgn = (nodesgn(end)+1):(nodesgn(end)+n(2)-1);
        gn((1:n(1):((n(2)-1)*n(1)))+n(1)-1,upperRight) = nodesgn';
        gn(((n(1)+1):n(1):(n(1)*n(2)))+n(1)-1,lowerRight) = nodesgn';
    end
    
    % domain corner nodes
    if nElements>1
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+4);
    else
        nodesgn = uint32(nInternalNodes+1):uint32(nInternalNodes+4);
    end  
    gn(1,lowerLeft) = nodesgn(1);
    gn(n(1),lowerRight) = nodesgn(2);
    gn((n(2)-1)*n(1)+1,upperLeft) = nodesgn(3);
    gn(n(1)*n(2),upperRight) = nodesgn(4);
  
    % remaining domain boundary nodes (the 'internal' ones)
    if p(1)>1
        % lower boundary
        nodesgn =  uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(p(1)-1));
        gn(1:n(1),bottomEdgeNodesInternalNumber) = reshape(nodesgn,[(p(1)-1) n(1)])';       
        % upper boundary
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(p(1)-1));
        gn(((n(2)-1)*n(1)+1):(n(1)*n(2)),topEdgeNodesInternalNumber) = reshape(nodesgn,[(p(1)-1) n(1)])';               
    end
    
    if p(2)>1
        % left boundary
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(p(2)-1));
        gn(1:n(1):((n(2)-1)*n(1))+1,leftEdgeNodesInternalNumber) = reshape(nodesgn,[(p(2)-1) n(2)])';        
        % right boundary
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(p(2)-1));
        gn((1:n(1):((n(2)-1)*n(1))+1)+n(1)-1,rightEdgeNodesInternalNumber) = reshape(nodesgn,[(p(2)-1) n(2)])';             
    end
    
end