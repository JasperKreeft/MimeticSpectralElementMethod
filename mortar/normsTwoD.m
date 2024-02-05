function normOutput = normsTwoD(mySol,funcPar,form,mapMethod,mapVar,errorQuad)
%
%   normsTwoD computes the error norm using gaussian quadrature in 2D
%
%   normOutput = normsTwoD(mySol,funcPar,form,mapMethod,mapVar,errorQuad)
%
%   input:
%       mySol       :: current solution
%       funcPar     :: function parameters 
%       form        :: 0 or 1 forms
%       mapMethod   :: 'scale' or 'deform'
%       mapVar      :: xBounds/yBounds for 'scale', mapPar for 'deform'
%       errorQuad   :: integration order (optional)
%
%   output:
%       normOutput  :: norm
%
%   Copyright 2011 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 09/11/2011 $     

%-------------------------------------------------------------------------%
% calculate norm                                                          %
%-------------------------------------------------------------------------%

    %%% 0-forms
    if form==0 

        if length(errorQuad)>2
            fprintf(':: errorQuad :: has too many values \n');
            return
        elseif length(errorQuad)==1
            errorQuad = [errorQuad errorQuad];
        end    
        if prod(errorQuad)==0
            fprintf(':: errorQuad :: has at least one zero \n');
            return
        end  
        
        % parameters                                                              
        px = size(mySol,2)-1;
        py = size(mySol,1)-1; 
        [quadNodesX quadWeightsX] = GaussQuad(errorQuad(1));     
        [quadNodesY quadWeightsY] = GaussQuad(errorQuad(2));           

        % interpolation                                             
        hnxQuad = MimeticpolyVal(quadNodesX,px,1);   
        hnyQuad = MimeticpolyVal(quadNodesY,py,1);  

        % mapping
        [quadNodesMapX quadNodesMapY] = mappingMethod(mapMethod,mapVar,quadNodesX,quadNodesY);  
        
        if strcmp('scale',mapMethod)
            J_GLL = jacobianTwoD(mapMethod,mapVar);
        elseif strcmp('deform',mapMethod)
            % J will not be a constant
            % NOTE: pay special attention to the multiplication by J (mat/vec)
        end
        J_Gauss = jacobianTwoD(mapMethod,mapVar); 
        % NOTE: for 'scale' J_GLL = J_Gauss = constant         

        interpSolXY = hnyQuad'*mySol.*J_GLL*hnxQuad;                           % interpolated solution
        exactSol = funcTwoD(quadNodesMapX,quadNodesMapY,0,funcPar,1).*J_Gauss; % exact solution
        diffSquared = (exactSol-interpSolXY).^2;                               % difference squared

        intX = zeros(errorQuad(1)+1,1);                            
        for i=1:errorQuad+1
            intX(i) = diffSquared(i,:)*quadWeightsX;	% apply weights in x-dir.
        end  
        intXY = intX'*quadWeightsY;	% apply weights in y-dir.
        
        normOutput = sqrt(intXY);	% norm output
 
    %%% 1-forms
    elseif form==1 
        % ...
    %%% 2-forms
    elseif form==2    
        normOutput = sqrt((quadTwoD(funcPar,0,mapVar,errorQuad)-sum(sum(mySol)))^2);
    end

end