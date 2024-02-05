function x = scaleMapOneD(xBounds,xi)
%
%   scaleMapOneD scales mesh nodes in 1D according to: 
%       x(xi)=a*xi+b (NOTE: affine transformation)
%
%   x = scaleMapOneD(xBounds,xi)
%
%   input:
%       xBounds :: = [xLeft xRight]
%       xi      :: reference coordinates
%
%   output:
%       x       :: mapped coordinates
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 01/10/2011 $  

%-------------------------------------------------------------------------%
% calculate mapping parameters                                            %
%-------------------------------------------------------------------------%

    a = (xBounds(2)-xBounds(1))/2; 
    b = xBounds(1)-(-1*a); 

%-------------------------------------------------------------------------%
% calculate mapped coordinates                                            %
%-------------------------------------------------------------------------%    
    
    x = a*xi+b;
    
end


 