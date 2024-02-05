function D = covariantDZero(p)
%dZero computes the discrete exterior derivative of zero forms.
%
%   D = dZero(p)
%
%   where D is the matrix of 1 and -1 entries that represents the discrete
%   exterior derivate (xi and eta components) of discrete zero forms. It is a
%   sparse matrix.
%
%   See also: DZEROXI, DZEROETA

%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 2010/01/21 13:49:00 $

    %% XI FINITE VOLUMES

    pX = p;
    pY = p+1;
    
    % xi part
    iiXiMinus = 1:(pX*(pY+1));
    jjXiMinus = 1:(pX*(pY+1));
    iiXiPlus = iiXiMinus;
    jjXiPlus = (1+pY+1):((pX+1)*(pY+1));
    
    % eta part
    iiEtaMinus = (1+pX*(pY+1)):(pX*(pY+1)+pY*(pX+1));
    jjEtaMinus = zeros([pY*(pX+1) 1]);
    
    for l = 1:(pX+1)
        for k = 1:pY
            jjEtaMinus((l-1)*pY + k) = (l-1)*(pY+1) + k;
        end
    end
    
    iiEtaPlus = iiEtaMinus;
    jjEtaPlus = jjEtaMinus + 1;
    
    % values
    valuesMinusXi = -ones([pX*(pY+1) 1]);
    valuesPlusXi = ones([pX*(pY+1) 1]);
    valuesMinusEta = -ones([pY*(pX+1) 1]);
    valuesPlusEta = ones([pY*(pX+1) 1]);
    
    D.Xi = sparse([iiXiMinus'; iiXiPlus'; iiEtaMinus'; iiEtaPlus'], [jjXiMinus'; jjXiPlus'; jjEtaMinus; jjEtaPlus], [valuesMinusXi; valuesPlusXi; valuesMinusEta; valuesPlusEta]);
    
    %% ETA FINITE VOLUMES

    pX = p+1;
    pY = p;
    
    % xi part
    iiXiMinus = 1:(pX*(pY+1));
    jjXiMinus = 1:(pX*(pY+1));
    iiXiPlus = iiXiMinus;
    jjXiPlus = (1+pY+1):((pX+1)*(pY+1));
    
    % eta part
    iiEtaMinus = (1+pX*(pY+1)):(pX*(pY+1)+pY*(pX+1));
    jjEtaMinus = zeros([pY*(pX+1) 1]);
    
    for l = 1:(pX+1)
        for k = 1:pY
            jjEtaMinus((l-1)*pY + k) = (l-1)*(pY+1) + k;
        end
    end
    
    iiEtaPlus = iiEtaMinus;
    jjEtaPlus = jjEtaMinus + 1;
    
    % values
    valuesMinusXi = -ones([pX*(pY+1) 1]);
    valuesPlusXi = ones([pX*(pY+1) 1]);
    valuesMinusEta = -ones([pY*(pX+1) 1]);
    valuesPlusEta = ones([pY*(pX+1) 1]);
    
    D.Eta = sparse([iiXiMinus'; iiXiPlus'; iiEtaMinus'; iiEtaPlus'], [jjXiMinus'; jjXiPlus'; jjEtaMinus; jjEtaPlus], [valuesMinusXi; valuesPlusXi; valuesMinusEta; valuesPlusEta]);
    
end