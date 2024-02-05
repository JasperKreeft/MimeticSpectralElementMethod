function gn = GlobalNumberingTwoFormTwoD(n,p)
%
%   GlobalNumberingTwoFormTwoD computes the global numbering of 2-form
%   degrees of freedom on a  primal mesh.
%
%   gn = GlobalNumberingTwoFormTwoD(n,p)
%
%   input:
%       n   :: the number of elements; if n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements
%       p   :: the order of the 0-form approximation; if p is a number then 
%              px=py=p, if it is a vector, then px!=py
%
%   output:
%       gn  :: global numbering of the surfaces of the mesh (gn is a N x p^2 
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
%       4.7 - bottom/top domain corner nodes
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
% parameters                                                              %
%-------------------------------------------------------------------------%   
    
    % compute the number of elements
    nElements = n(1)*n(2);

    % compute the total number of degrees of freedom
    dof = p(1)*p(2)*nElements;

%-------------------------------------------------------------------------%
% compute global numbering                                                %
%-------------------------------------------------------------------------%    
    
    gn = reshape(uint32(1):uint32(dof),p(1)*p(2),nElements)';

end