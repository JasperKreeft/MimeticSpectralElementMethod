function [J D] = jacobianTwoD(mapMethod,mapVar,xi,eta)
%
%   jacobianTwoD is the general jacobian function in 2D
%
%   [J D] = jacobianTwoD(mapMethod,mapVar,xi,eta)
%
%   input:
%       mapMethod   :: 'scale' or 'deform'
%       mapVar      :: xBounds/yBounds for 'scale', mapPar for 'deform'
%       xi          :: reference coordinates in x-dir.
%       eta         :: reference coordinates in y-dir.
%
%   output:
%       J           :: jacobian
%       D           :: derivatives
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 25/10/2011 $  

%-------------------------------------------------------------------------%
% input correction                                                        %
%-------------------------------------------------------------------------%

    if nargin>2
        if size(xi,1)<size(xi,2)
            xi = xi';
        end
        if size(eta,1)<size(eta,2)
            eta = eta';
        end
    end

%-------------------------------------------------------------------------%
% scale                                                                   %
%-------------------------------------------------------------------------%

    if strcmp('scale',mapMethod)
        if nargout>1
            [J D] = scaleMapJacTwoD(mapVar);
        else
            J = scaleMapJacTwoD(mapVar);
        end
    end

%-------------------------------------------------------------------------%
% deform                                                                  %
%-------------------------------------------------------------------------%

    if strcmp('deform',mapMethod)
        fprintf(':: mapMethod :: map method >deform< not implemented \n')  
    end      

end