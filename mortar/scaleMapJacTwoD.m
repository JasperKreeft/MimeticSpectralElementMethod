function [J D] = scaleMapJacTwoD(mapVar)
%
%   scaleMapJacTwoD calculates the 2D jacobian of the prescribed scaling: 
%       x(xi)=a*xi+b    (NOTE: affine transformation)
%       y(eta)=c*eta+d  (NOTE: affine transformation)
%
%   [J D] = scaleMapJacTwoD(mapVar)
%
%   input:
%       mapVar  :: xBounds/yBounds for 'scale', mapPar for 'deform'
%
%   output:
%       J    	:: jacobian
%       D      	:: derivatives
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 25/10/2011 $  

%-------------------------------------------------------------------------%
% bounds                                                                  %
%-------------------------------------------------------------------------%

    xBounds = mapVar(1:2);
    yBounds = mapVar(3:4);

%-------------------------------------------------------------------------%
% mapping parameters                                                      %
%-------------------------------------------------------------------------%

    a = (xBounds(2)-xBounds(1))/2; 
    % b = xBounds(1)-(-1*a); 
    c = (yBounds(2)-yBounds(1))/2;
    % d = xBounds(1)-(-1*a);     

%-------------------------------------------------------------------------%
% calculate jacobian                                                      %
%-------------------------------------------------------------------------%    
    
    dx_dxi = a;             % derivative    
    dy_deta = c;            % derivative 
    J = dx_dxi*dy_deta;     % calculate jacobian
    if nargout>1
        D = zeros(2);       % matrix of derivatives
        D(1,1) = dx_dxi;    % fill matrix
        D(2,2) = dy_deta;   % fill matrix
    end

    
end





 