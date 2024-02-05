function D = dZeroEta(p)
%dZeroXi computes the Eta component of the discrete exterior derivative of zero forms.
%
%   D = dZeroEta(p)
%
%   where D is the matrix of 1 and -1 entries that represents the discrete
%   exterior derivate (eta component) of discrete zero forms. It is a
%   sparse matrix.
%
%   See also: DZERO, DZEROXI

%   Copyright 2010 Artur Palha
%   $Revision: 1.0 $  $Date: 2010/01/21 13:49:00 $

    % eta part
    iiEtaMinus = 1:(p*(p+1));
    jjEtaMinus = zeros([p*(p+1) 1]);
    
    for l = 1:(p+1)
        for k = 1:p
            jjEtaMinus((l-1)*p + k) = (l-1)*(p+1) + k;
        end
    end
    
    iiEtaPlus = iiEtaMinus;
    jjEtaPlus = jjEtaMinus + 1;
    
    % values
    valuesMinus = -ones([p*(p+1) 1]);
    valuesPlus = ones([p*(p+1) 1]);
    
    D = sparse([iiEtaMinus'; iiEtaPlus'], [jjEtaMinus; jjEtaPlus], [valuesMinus; valuesPlus]);
end