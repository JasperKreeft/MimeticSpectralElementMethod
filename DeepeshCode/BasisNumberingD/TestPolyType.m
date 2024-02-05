function result = TestPolyType(polyType)
%TestPolyType Returns true if the polyType is a valid one.
%   result = TestPolyType(polyType) checks if the polyType is one of the
%   valid types of polynomials. Valid types:
%
%       Lobatto, Gauss, EGauss
%
%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/11/27 13:38:00 $

    validPolyTypes = strvcat('Lobatto', 'Gauss', 'EGauss');
    
    sizePolyType = size(polyType);
    
    result = true;
    
    for n=1:sizePolyType(1)
        if isempty(strmatch(polyType(n,:), validPolyTypes,'exact')), 
            result = result && false;
        else
            result = result && true;
        end
    end
end