function D = dZeroXiTwoD(p)
%
%   dZeroXi computes the xi component of the discrete exterior derivative 
%   of zero forms (local (per element) numbering applies)
%
%   D = dZeroXiTwoD(p)
%
%   input:
%       p   :: the order of the 0-form approximation; if p is a number then 
%              px=py=p, if it is a vector, then px!=py
%
%   output:
%       D   :: the matrix of 1 and -1 entries that represents the discrete
%              exterior derivate (xi part) of discrete zero forms (sparse)
%
%   Copyright 2010 Artur Palha
%   $Revision: 1.0 $  $Date: 21/01/2010 $
%   Revised by: Peter Kuystermans
%   $Revision: 2.0 $  $Date: 13/07/2011 $  
%   Revised by: Peter Kuystermans
%   $Revision: 3.0 $  $Date: 07/11/2011 $ 

%-------------------------------------------------------------------------%
% input checks                                                            %
%-------------------------------------------------------------------------%   

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
% incidence matrix                                                        %
%-------------------------------------------------------------------------%   

    % xi part
    iiXiMinus = (1:(p(1)*(p(2)+1)))'; % line (-)
    iiXiPlus = iiXiMinus; % line (+)
    jjXiMinus = zeros([p(1)*(p(2)+1) 1]); % point (-)
    for i=1:p(2)+1
        jjXiMinus((1:p(1))+p(1)*(i-1)) = (1:p(1))+(p(1)+1)*(i-1);
    end
    jjXiPlus = jjXiMinus+1; % point (+)
 
    % values
    valuesMinus = -ones([p(1)*(p(2)+1) 1]);
    valuesPlus = ones([p(1)*(p(2)+1) 1]);
    
    % building D
    D = sparse([iiXiMinus;iiXiPlus],[jjXiMinus;jjXiPlus],[valuesMinus;valuesPlus]);
    
end