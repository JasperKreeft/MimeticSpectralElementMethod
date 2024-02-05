function D = covariantDOne(p, varargin)
%dOne computes the discrete exterior derivative of one forms.
%
%   D = dOne(p)
%
%   where D is the matrix of 1 and -1 entries that represents the discrete
%   exterior derivate (xi and eta part) of discrete one forms. It is a
%   sparse matrix.
%
%   See also: DONEXI, DONEETA

%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 2012/4/17$

    if (size(varargin,2))
        connection = varargin{1};
    else
        connection = 0;
    end
    
    %% D for Finite volumes surrounding Xi Edges
    
    pX = p;
    pY = p+1;
    
    % eta part
    iiEtaMinus = (1:(pX*pY))';
    jjEtaMinus = ((1:(pX*pY)) + pX*(pY+1))';
    iiEtaPlus = iiEtaMinus;
    jjEtaPlus = jjEtaMinus + pY;
    
    % xi part
    iiXiMinus = iiEtaMinus;
    iiXiPlus = iiXiMinus;
    jjXiPlus = repmat((1:pY)',1,pX) + repmat(0:(pY+1):(pX-1)*(pY+1),pY,1);
    jjXiPlus = jjXiPlus(:);
    jjXiMinus = jjXiPlus + 1;
    
    % values
    valuesMinus = -ones([pX*pY 1]);
    valuesPlus = ones([pX*pY 1]);
    
    D.Xi = sparse([iiXiMinus; iiXiPlus; iiEtaMinus; iiEtaPlus], [jjXiMinus; jjXiPlus; jjEtaMinus; jjEtaPlus], [valuesMinus; valuesPlus; valuesMinus; valuesPlus], pX*pY, pX*(pY+1)+pY*(pX+1));
    
    %% D for Finite volumes surrounding Eta Edges
    
    pX = p+1;
    pY = p;
    
    % eta part
    iiEtaMinus = (1:(pX*pY))';
    jjEtaMinus = ((1:(pX*pY)) + pX*(pY+1))';
    iiEtaPlus = iiEtaMinus;
    jjEtaPlus = jjEtaMinus + pY;
    
    % xi part
    iiXiMinus = iiEtaMinus;
    iiXiPlus = iiXiMinus;
    jjXiPlus = repmat((1:pY)',1,pX) + repmat(0:(pY+1):(pX-1)*(pY+1),pY,1);
    jjXiPlus = jjXiPlus(:);
    jjXiMinus = jjXiPlus + 1;
    
    % values
    valuesMinus = -ones([pX*pY 1]);
    valuesPlus = ones([pX*pY 1]);
    
    D.Eta = sparse([iiXiMinus; iiXiPlus; iiEtaMinus; iiEtaPlus], [jjXiMinus; jjXiPlus; jjEtaMinus; jjEtaPlus], [valuesMinus; valuesPlus; valuesMinus; valuesPlus], pX*pY, pX*(pY+1)+pY*(pX+1));
    
    %% Implement Connection Forms
    
    
    
end