function innerProduct = OneFormInnerOneForm1D(p,pint,g,gridType,element)
% ZeroFormInnerZeroForm1D computes the inner product between the
% 1-form basis functions and 0-form basis functions of order p, with
% integration order pInt and metric g.
%
% innerProduct = OneFormInnerOneForm1D(p,pInt,g,gridType)
%
%   INPUTS
%
%       p       :: the order of the basis functions
%       pint    :: the order of integration
%       g       :: the metric evaluated at the integration points
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 2011/02/21$

    % compute the integration nodes and the integration weights
    [quadNodes quadWeights] = GaussQuad(pint);
        
    % compute the 1-form basis functions, evaluated at the integrations
    % nodes
    oneFormBasis = EdgeFunction(quadNodes,p,gridType);
        
    % Evaluate metric
    gEvaluated = g{element}(quadNodes);

    % compute the inner product between the 1-form basis functions
    innerProduct = oneFormBasis*spdiags(quadWeights(:)./gEvaluated(:),0,pint+1,pint+1)*oneFormBasis';

end