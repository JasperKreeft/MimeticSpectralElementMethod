function D = dOne(p)
%dOne computes the discrete exterior derivative of one forms.
%
%   D = dOne(p)
%
%   where D is the matrix of 1 and -1 entries that represents the discrete
%   exterior derivate (xi and eta part) of discrete one forms. It is a
%   sparse matrix.
%
%   See also: DONEXI, DONEETA

%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/12/17 13:49:00 $
    
    % eta part
    iiEtaMinus = (1:(p*p))';
    jjEtaMinus = ((1:(p*p)) + p*(p+1))';
    iiEtaPlus = iiEtaMinus;
    jjEtaPlus = jjEtaMinus + p;
    
    % xi part
    iiXiMinus = iiEtaMinus;
    iiXiPlus = iiXiMinus;
    jjXiPlus = zeros([p*p 1]);
    
    for l = 1:p
        for k = 1:p
            jjXiPlus((l-1)*p + k) = (l-1)*(p+1) + k;
        end
    end
    
    jjXiMinus = jjXiPlus + 1;
    
    % values
    valuesMinus = -ones([p*p 1]);
    valuesPlus = ones([p*p 1]);
    
    D = sparse([iiXiMinus; iiXiPlus; iiEtaMinus; iiEtaPlus], [jjXiMinus; jjXiPlus; jjEtaMinus; jjEtaPlus], [valuesMinus; valuesPlus; valuesMinus; valuesPlus], p^2, 2*p*(p+1));
end