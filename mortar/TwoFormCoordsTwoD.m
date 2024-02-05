function [X Y] = TwoFormCoordsTwoD(n,p,x,y)
%
%   TwoFormCoordsTwoD(n,p,x,y) generates the two form coordinates of a 2D mesh 
%   consisting of n elements and order p
%
%   [X Y] = TwoFormCoordsTwoD(n,p,x,y)
%
%   input:
%       n       :: is the number of elements (if n is a number then n x n elements
%                  are assumed, if n is a vector with two entries then n(1) x n(2)
%                  elements are assumed)
%       p       :: the order of the 2-form approximation; if p is a number then 
%                  px=py=p, if it is a vector, then px!=py
%       x       :: contains the x boundaries of the mesh
%       y       :: contains the y boundaries of the mesh
%
%   output:
%       X       :: x-coordinates (X is an N x [px*py] x 2 matrix where N is the total 
%                  number of elements and [px*py] the number of 2-forms per element; 
%                  the x-coordinates of the corner nodes LL, LR, UL, UR (i.e. local numbering)  
%                  are stored at levels X(:,:,1), X(:,:,2), X(:,:,3) and X(:,:,4) respectively 
%       Y       :: y-coordinates (Y is an N x [px*py] x 2 matrix where N is the total 
%                  number of elements and [px*py] the number of 2-forms per element; 
%                  the y-coordinates of the corner nodes LL, LR, UL, UR (i.e. local numbering)  
%                  are stored at levels Y(:,:,1), Y(:,:,2), Y(:,:,3) and Y(:,:,4) respectively 
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 21/10/2010 $

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
% mesh parameters                                                         %
%-------------------------------------------------------------------------%     
    
    % compute the number of elements
    nElements = n(1)*n(2);
    nElementSurfaces = p(1)*p(2);
    % nodes per element
%     nElementNodes = (p(1)+1)*(p(2)+1);

%-------------------------------------------------------------------------%
% storage                                                                 %
%-------------------------------------------------------------------------%    
       
    lowerLeftNodes = zeros(n(1)*n(2),2);        % lower left element nodes    
    X = zeros(nElements,nElementSurfaces,4);    % x-coordinates
    Y = zeros(nElements,nElementSurfaces,4);    % y-coordinates    

%-------------------------------------------------------------------------%
% calculated parameters                                                   %
%-------------------------------------------------------------------------%
    
    % compute domain size
    Lx = x(2)-x(1);
    Ly = y(2)-y(1);
    
    % element lengths
    deltax = Lx/n(1); 
    deltay = Ly/n(2); 
    
    % compute the x coordinates of the lower left element nodes
    lowerLeftNodes(:,1) = repmat(x(1)+(0:(n(1)-1))*deltax,[1 n(2)])';
    
    % compute the y coordinates of the lower left nodes
    lowerLeftNodes(:,2) = rectpulse(y(1)+(0:(n(2)-1))*deltay,n(1))';  
    
    % generate element nodes 
    xi = GLLnodes(p(1));
    eta = GLLnodes(p(2));

%-------------------------------------------------------------------------%
% matrix build-up                                                         %
%-------------------------------------------------------------------------%
    
    % loop over the elements and plot all the nodes and edges
    for k=1:nElements
        
        % scale the nodes
        x = (0.5*(xi+1)*deltax+lowerLeftNodes(k,1))';
        y = (0.5*(eta+1)*deltay+lowerLeftNodes(k,2))';
        
        % x-coordinates
        X(k,:,1) = repmat(x(1:end-1),1,p(2));
        X(k,:,2) = repmat(x(2:end),1,p(2));
        X(k,:,3:4) = X(k,:,1:2);

        % y-coordinates
        Y(k,:,1) = reshape(repmat(y(1:end-1),p(1),1),1,nElementSurfaces);
        Y(k,:,2) = Y(k,:,1);
        Y(k,:,3) = reshape(repmat(y(2:end),p(1),1),1,nElementSurfaces); 
        Y(k,:,4) = Y(k,:,3);
        
    end
    
end