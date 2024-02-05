function [x y] = mappingMethod(mapMethod,mapVar,xi,eta)
%
%   mappingOneD is the general transformation function in 1D and 2D
%
%   [x y] = mappingMethod(mapMethod,mapVar,xi,eta)
%
%   input:
%       mapMethod   :: 'scale' or 'deform'
%       mapVar      :: xBounds/yBounds for 'scale', mapPar for 'deform'
%       xi          :: reference coordinates in x-dir.
%       eta         :: reference coordinates in y-dir.
%
%   output:
%       x           :: mapped coordinates in x-dir.
%       y           :: mapped coordinates in y-dir.
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 25/10/2011 $  

%-------------------------------------------------------------------------%
% input correction                                                        %
%-------------------------------------------------------------------------%

    if size(xi,1)<size(xi,2)
        xi = xi';
    end
    if nargin>3
        if size(eta,1)<size(eta,2)
            eta = eta';
        end
    end

%-------------------------------------------------------------------------%
% scale                                                                   %
%-------------------------------------------------------------------------%

    if strcmp('scale',mapMethod)
        xBounds = mapVar(1:2);
        x = scaleMapOneD(xBounds,xi);
        if nargout>1
            if length(mapVar)~=4
                fprintf(':: mapVar :: expected length is 4, got 2 \n')
                fprintf(':: mapVar :: using mapVar(1:2) \n')  
                mapVarRep = mapVar;
                mapVar = zeros(4,1);
                mapVar(1:2) = mapVarRep;
                mapVar(3:4) = mapVarRep;
            else
                if nargin<4
                    fprintf(':: eta :: additional reference coordinates expected \n')
                    fprintf(':: eta :: using xi \n')
                    eta = xi;                          
                end
            end
            yBounds = mapVar(3:4);
            y = scaleMapOneD(yBounds,eta);  
        end
    end

%-------------------------------------------------------------------------%
% deform                                                                  %
%-------------------------------------------------------------------------%

    if strcmp('deform',mapMethod)
        fprintf(':: mapMethod :: map method >deform< not implemented \n')        
    end    
        
end