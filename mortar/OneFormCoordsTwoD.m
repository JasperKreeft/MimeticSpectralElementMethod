function [X Y] = OneFormCoordsTwoD(n,p,x,y)
%
%   OneFormCoordsTwoD generates the zero form coordinates of a 2D mesh 
%   consisting of n elements and order p
%
%   [X Y] = OneFormPrimal(n,p,x,y)
%
%   input:
%       n       :: is the number of elements (if n is a number then n x n elements
%                  are assumed, if n is a vector with two entries then n(1) x n(2)
%                  elements are assumed)
%       p       :: the order of the 0-form approximation; if p is a number then 
%                  px=py=p, if it is a vector, then px!=py
%       x       :: contains the x boundaries of the mesh
%       y       :: contains the y boundaries of the mesh
%
%   output:
%       X       :: x-coordinates (X is an N x nEdges x 2 matrix where N is the total 
%                  number of elements and nEdges the number of edges [px*(py+1)+(px+1)*py]; 
%                  the left/lower edge x-coordinate is stored at level X(:,:,1) and the 
%                  right/upper edge x-coordinate at level X(:,:,2)) 
%       Y       :: y-coordinates (Y is an N x nEdges x 2 matrix where N is the total 
%                  number of elements and nEdges the number of edges [px*(py+1)+(px+1)*py]; 
%                  the left/lower edge y-coordinate is stored at level Y(:,:,1) and the 
%                  right/upper edge y-coordinate at level Y(:,:,2)) 
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
    
    % compute the number of edges in each element
    nElementEdgesX = p(1)*(p(2)+1);
    nElementEdgesY = p(2)*(p(1)+1);
    
%-------------------------------------------------------------------------%
% storage                                                                 %
%-------------------------------------------------------------------------%  

    lowerLeftNodes = zeros(n(1)*n(2),2);        % lower left elements nodes
    Xxi = zeros(nElements,nElementEdgesX,2);    % x-coordinates h-edges  
    Yxi = zeros(nElements,nElementEdgesX,2);    % y-coordinates h-edges
    Xeta = zeros(nElements,nElementEdgesY,2);	% x-coordinates v-edges  
    Yeta = zeros(nElements,nElementEdgesY,2);   % y-coordinates v-edges    
    
%-------------------------------------------------------------------------%
% calculated parameters                                                   %
%-------------------------------------------------------------------------%   
   
    % compute domain size
    Lx = x(2)-x(1);
    Ly = y(2)-y(1);
    
    % element lengths
    deltax = Lx/n(1); 
    deltay = Ly/n(2); 
    
    % compute the x coordinates of the lower left nodes
    lowerLeftNodes(:,1) = repmat(x(1)+(0:(n(1)-1))*deltax,[1 n(2)])';
    
    % compute the y coordinates of the lower left nodes
    lowerLeftNodes(:,2) = rectpulse(y(1)+(0:(n(2)-1))*deltay,n(1))';
    
    % generate element nodes
    xi = GLLnodes(p(1));
    eta = GLLnodes(p(2));
    
%-------------------------------------------------------------------------%
% coordinate calculation and matrix build-up                              %
%-------------------------------------------------------------------------%

    for k=1:nElements
   
        % scale the nodes
        x = (0.5*(xi+1)*deltax+lowerLeftNodes(k,1))';
        y = (0.5*(eta+1)*deltay+lowerLeftNodes(k,2))';
               
        % horizontal edges
        Xxi(k,:,1) = repmat(x(1:end-1),1,p(2)+1);
        Xxi(k,:,2) = repmat(x(2:end),1,p(2)+1);
        Yxi(k,:,1) = reshape(repmat(y,p(1),1),1,p(1)*(p(2)+1));
        Yxi(k,:,2) = reshape(repmat(y,p(1),1),1,p(1)*(p(2)+1));
        
        % vertical edges
        Xeta(k,:,1) = reshape(repmat(x,p(2),1),1,(p(1)+1)*p(2));
        Xeta(k,:,2) = reshape(repmat(x,p(2),1),1,(p(1)+1)*p(2));
        Yeta(k,:,1) = repmat(y(1:end-1),1,p(1)+1);
        Yeta(k,:,2) = repmat(y(2:end),1,p(1)+1);
       
    end 
    
    % matrix build-up
    X = [Xxi Xeta];
    Y = [Yxi Yeta];
    
end