function D = dZeroXi(p)
%dZeroXi computes the Xi component of the discrete exterior derivative of zero forms.
%
%   D = dZeroXi(p)
%
%   where D is the matrix of 1 and -1 entries that represents the discrete
%   exterior derivate (xi component) of discrete zero forms. It is a
%   sparse matrix.
%
%   See also: DZERO, DZEROETA

%   Copyright 2010 Artur Palha
%   $Revision: 1.0 $  $Date: 2010/01/21 13:49:00 $

    % xi part
    iiXiMinus = 1:(p*(p+1));
    jjXiMinus = 1:(p*(p+1));
    iiXiPlus = iiXiMinus;
    jjXiPlus = (1+p+1):((p+1)*(p+1));
    
    % values
    valuesMinus = -ones([p*(p+1) 1]);
    valuesPlus = ones([p*(p+1) 1]);
    
    D = sparse([iiXiMinus'; iiXiPlus'], [jjXiMinus'; jjXiPlus'], [valuesMinus; valuesPlus]);
end