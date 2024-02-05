function quadVal = quadTwoD(funcPar,fD,mapVarCell,pInt,mapVar)
%
%   quadTwoD calculates the integral of a function in 2D
%
%   quadTwoD(funcPar,fD,mapVarCell,pInt,mapVar)
%
%   input:
%       funcPar     :: function parameters
%       fD          :: derivative
%       mapVarCell  :: integration interval (cell)
%       pInt        :: integration order [pIntX pIntY]
%       mapVar      :: integration interval (domain)
%
%   output:
%       quadVal :: integral value
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 23/11/2011 $  

%-------------------------------------------------------------------------%
% input checks                                                            %
%-------------------------------------------------------------------------%

    if length(pInt)>2
        fprintf(':: p :: has too many values \n');
        return
    elseif length(pInt)==1
        pIntX = pInt;
        pIntY = pInt;
    elseif length(pInt)==2 
        pIntX = pInt(1);
        pIntY = pInt(2);
    end
    
    if prod(pInt)==0
        fprintf(':: p :: has at least one zero \n');
        return
    end  

%-------------------------------------------------------------------------%
% function evaluations provided                                           %
%-------------------------------------------------------------------------%
 
...
 
%-------------------------------------------------------------------------%
% function evaluations not provided                                       %
%-------------------------------------------------------------------------%       

    %%% quadrature                                                      
    [quadNodesX quadWeightsX] = GaussQuad(pIntX);
    [quadNodesY quadWeightsY] = GaussQuad(pIntY);

    %%% function evaluation
    if nargin < 5
        
        %%% mapping
        [quadNodesMappedX quadNodesMappedY] = ...
            mappingMethod('scale',mapVarCell,quadNodesX,quadNodesY);
        J = jacobianTwoD('scale',mapVarCell,quadNodesX,quadNodesY);        
        
        fVal = funcTwoD(quadNodesMappedX,quadNodesMappedY,fD,funcPar,1); 
   
    else
        
        %%% mapping
        [quadNodesMappedX quadNodesMappedY] = ...
            mappingMethod('scale',mapVarCell,quadNodesX,quadNodesY);
        J = jacobianTwoD('scale',mapVarCell,quadNodesX,quadNodesY);        
        
        fVal = particularSolutionTwoD(funcPar,fD,mapVar(1:2),mapVar(3:4),quadNodesMappedX,quadNodesMappedY); 
  
    end
    
    %%% apply quadrature
    intX = fVal*quadWeightsX;
    intY = quadWeightsY'*intX;
    intMap = intY*J;
    
    %%% output
    quadVal = intMap;

end