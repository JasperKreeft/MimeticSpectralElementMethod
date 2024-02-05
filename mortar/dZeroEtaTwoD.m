function D = dZeroEtaTwoD(p)
%
%   dZeroEta computes the eta component of the discrete exterior derivative 
%   of zero forms (local (per element) numbering applies)
%
%   D = dZeroEtaTwoD(p)
%
%   input:
%       p   :: the order of the 0-form approximation; if p is a number then 
%              px=py=p, if it is a vector, then px!=py
%
%   output:
%       D   :: the matrix of 1 and -1 entries that represents the discrete
%              exterior derivate (eta part) of discrete zero forms (sparse)
%
%   Copyright 2010 Artur Palha
%   $Revision: 1.0 $  $Date: 21/01/2010 $
%   Revised by: Peter Kuystermans
%   $Revision: 2.0 $  $Date: 13/07/2011 $   
%   Revised by: Peter Kuystermans
%   $Revision: 3.0 $  $Date: 18/10/2011 $ 

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

    % eta part
    iiEtaMinus = (1:(p(2)*(p(1)+1)))'; % line(-)  
    iiEtaPlus = iiEtaMinus; % line(+)    
    jjEtaMinus = zeros([p(2)*(p(1)+1) 1]); % point(-)
    for i=1:p(1)+1    
        jjEtaMinus((1:p(2))+p(2)*(i-1)) = (1:(p(1)+1):(1+(p(1)+1)*(p(2)-1)))+(i-1);
    end
    jjEtaPlus = jjEtaMinus+(p(1)+1); % point(+)

    % values
    valuesMinus = -ones([p(2)*(p(1)+1) 1]);
    valuesPlus = ones([p(2)*(p(1)+1) 1]);

    % building D
    D = sparse([iiEtaMinus;iiEtaPlus],[jjEtaMinus;jjEtaPlus],[valuesMinus;valuesPlus]);

end