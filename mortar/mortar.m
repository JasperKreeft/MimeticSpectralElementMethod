clear all
close all
clc

%%% ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo %%%
%%% mortar element test case for 0-forms in 2d                          %%%
%%% p_left = left element order, x_left/y_left = element dimensions     %%%
%%% p_right = right element order, x_right/y_right = element dimensions %%% 
%%% ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo %%%

%-------------------------------------------------------------------------%
% settings                                                                %
%-------------------------------------------------------------------------%

% element order
orderLeft = 2:3;  % left element order
orderRight = 2:3; % right element order

% number of elements
hLeft = 1;  % number of left elements (NOTE: not implemented yet)
hRight = 1; % number of right elements (NOTE: not implemented yet)

% global domain
mapMethod = 'scale';	% mapping method (NOTE: only use 'scale')
mapVar = [-1 1 -1 1];	% mapping variables; [xLeft xRight yLeft yRight]

% left element
mapMethodLeft = 'scale';	% mapping method (NOTE: only use 'scale')    
mapVarLeft = [-1 0 -1 1];   % mapping variables; [xLeft xRight yLeft yRight]	

% right element
mapMethodRight = 'scale';	% mapping method (NOTE: only use 'scale') 
mapVarRight = [0 1 -1 1];   % mapping variables; [xLeft xRight yLeft yRight] 	

funcPar = [3 1 1];	% [function, function parameters wx/wy]
fineGridN = 100;   	% refinement level of fine grid
errorQuad = 100;   	% integration order norm quadrature

BCs = [1 1 1 1];	% [S(1) N(2) W(3) E(4)], 1 = dirichlet, 2 = neumann

%-------------------------------------------------------------------------%
% variable initialization                                                 %
%-------------------------------------------------------------------------%

figureCount = 1;	% figure counter

%-------------------------------------------------------------------------%
% bc checks                                                               %
%-------------------------------------------------------------------------%

dirichletSides = find(BCs);	% sides with dirichlet bc
neumannSides = find(BCs-1); % sides with neumann bc
dirichletFlag = 1;          % dirichlet flag
if isempty(dirichletSides)  % if empty then...
	dirichletFlag = 0;      % ...no dirichlet bc
end
neumannFlag = 1;            % dirichletflag
if isempty(neumannSides)    % if empty then...
    neumannFlag = 0;        % ...no neumann bc
end
    
%-------------------------------------------------------------------------%
% domain definition                                                       %
%-------------------------------------------------------------------------%

%%% left element dimensions
xElementCoordsLeft = mapVarLeft(1:2);   
yElementCoordsLeft = mapVarLeft(3:4);  

%%% right element dimensions
xElementCoordsRight = mapVarRight(1:2);   
yElementCoordsRight = mapVarRight(3:4);

% check element connections (dimension wise)
if xElementCoordsLeft(2)~=xElementCoordsRight(1)
    fprintf('Check input: Elements not connected in x-direction \n');
	return
end
if (yElementCoordsLeft(1)~=yElementCoordsRight(1))||(yElementCoordsLeft(2)~=yElementCoordsRight(2))
    fprintf('Check input: Elements not connected in y-direction \n');
	return
end

% % 0-forms
% figure(figureCount)
% % subplot(1,2,[1 2])
% subplot(1,2,1), hold on
% PlotGlobalNumberingZeroFormTwoD(hLeft,pLeft,xElementCoordsLeft,yElementCoordsLeft,1,1)
% subplot(1,2,2), hold on
% PlotGlobalNumberingZeroFormTwoD(hRight,pRight,xElementCoordsRight,yElementCoordsRight,1,1)
% figureCount = figureCount+1;

%-------------------------------------------------------------------------%
% calculate exact solution (fine grid)                                    %
%-------------------------------------------------------------------------%

fineGridRef = linspace(-1,1,fineGridN);                                     % finegrid    
[xFine yFine] = mappingMethod(mapMethod,mapVar,fineGridRef,fineGridRef);	% interpolation mesh x-dir./y-dir.
U_Fine = funcTwoD(xFine,yFine,0,funcPar,1);                                 % exact solution (matrix) 
F_Fine = -funcTwoD(xFine,yFine,2,funcPar,1);

for pLeft=orderLeft
    for pRight=orderRight

        %-------------------------------------------------------------------------%
        % standard element mesh and polynomials                                   %
        %-------------------------------------------------------------------------%

        %%% left element
        [refMeshNodesLeft refMeshWeightsLeft] = GLLnodes(pLeft);	% mesh reference nodes and weights (x-dir.)
        [hnLeft enLeft] = MimeticpolyVal(refMeshNodesLeft,pLeft,1);	% coarse mesh interpolants (x-dir.) 
        [hfLeft efLeft] = MimeticpolyVal(fineGridRef,pLeft,1);      % fine mesh interpolants (x-dir.)

        %%% right element
        [refMeshNodesRight refMeshWeightsRight] = GLLnodes(pRight);     % mesh reference nodes and weights (x-dir.)
        [hnRight enRight] = MimeticpolyVal(refMeshNodesRight,pRight,1);	% coarse mesh interpolants (x-dir.) 
        [hfRight efRight] = MimeticpolyVal(fineGridRef,pRight,1);       % fine mesh interpolants (x-dir.)

        %-------------------------------------------------------------------------%
        % mesh topology                                                           % 
        %-------------------------------------------------------------------------%    

        %%% left element
        D10xLeft = dZeroXiTwoD(pLeft);	% incidence matrix x-dir. 
        D10yLeft = dZeroEtaTwoD(pLeft);	% incidence matrix y-dir.

        ZeroFormsLeft = GlobalNumberingZeroFormTwoD(hLeft,pLeft);	% 0-forms local numbering
        OneFormsLeft  = GlobalNumberingOneFormTwoD(hLeft,pLeft);	% 1-forms local numbering
        TwoFormsLeft  = GlobalNumberingTwoFormTwoD(hLeft,pLeft);  	% 2-forms local numbering    

        nNodesLeft = double(max(max(ZeroFormsLeft)));	% # of 0-forms
        nLinesLeft = double(max(max(OneFormsLeft)));   	% # of 1-forms
        nSurfacesLeft = double(max(max(TwoFormsLeft)));	% # of 2-forms

        %%% right element
        D10xRight = dZeroXiTwoD(pRight);	% incidence matrix x-dir.
        D10yRight = dZeroEtaTwoD(pRight);	% incidence matrix y-dir.

        ZeroFormsRight = GlobalNumberingZeroFormTwoD(hRight,pRight);	% 0-forms local numbering
        OneFormsRight  = GlobalNumberingOneFormTwoD(hRight,pRight);     % 1-forms local numbering
        TwoFormsRight  = GlobalNumberingTwoFormTwoD(hRight,pRight);  	% 2-forms local numbering    

        nNodesRight = double(max(max(ZeroFormsRight)));     % # of 0-forms
        nLinesRight = double(max(max(OneFormsRight)));   	% # of 1-forms
        nSurfacesRight = double(max(max(TwoFormsRight)));	% # of 2-forms

        % total numbers
        nNodes = nNodesLeft+nNodesRight;
        nLines = nLinesLeft+nLinesRight;
        nSurfaces = nSurfacesLeft+nSurfacesRight;

        %-------------------------------------------------------------------------%
        % element submatrices                                                     %
        %-------------------------------------------------------------------------%

        %%% left element
        % mesh
        [xMeshLeft yMeshLeft]= mappingMethod(mapMethodLeft,mapVarLeft,refMeshNodesLeft,refMeshNodesLeft);
        % RHS
        fMeshLeft = -funcTwoD(xMeshLeft,yMeshLeft,2,funcPar,2);
        % jacobian                           
        [Jleft Dleft] = jacobianTwoD(mapMethodLeft,mapVarLeft);	% NOTE: J is just a constant for 'scale'
        dxdxiLeft = Dleft(1,1);                           
        dydetaLeft = Dleft(2,2);
        % matrix construction
        Wleft = spdiags(refMeshWeightsLeft,0,pLeft+1,pLeft+1);                  % integration weights x-dir.
        M0Left = hnLeft*Wleft*hnLeft';                                          % lagrange integration x-dir.
        M1Left = enLeft*Wleft*enLeft';                                          % edge integration x-dir.
        A1Left = D10xLeft'*kron(M0Left,M1Left)*D10xLeft/Jleft*dydetaLeft^2;	% contribution in x-dir.
        A2Left = D10yLeft'*kron(M0Left,M1Left)*D10yLeft/Jleft*dxdxiLeft^2;     % contribution in y-dir.
        RHSleft = kron(M0Left,M0Left)*fMeshLeft*Jleft;                          % right hand side
        Aleft = (A1Left+A2Left);

        %%% right element
        % mesh
        [xMeshRight yMeshRight] = mappingMethod(mapMethodRight,mapVarRight,refMeshNodesRight,refMeshNodesRight);
        % RHS
        fMeshRight = -funcTwoD(xMeshRight,yMeshRight,2,funcPar,2);
        % jacobian                           
        [Jright Dright] = jacobianTwoD(mapMethodRight,mapVarRight);	% NOTE: J is just a constant for 'scale'
        dxdxiRight = Dright(1,1);                           
        dydetaRight = Dright(2,2);  
        % matrix construction
        Wright = spdiags(refMeshWeightsRight,0,pRight+1,pRight+1);                  % integration weights x-dir.
        M0Right = hnRight*Wright*hnRight';                                          % lagrange integration x-dir.
        M1Right = enRight*Wright*enRight';                                          % edge integration x-dir.
        A1Right = D10xRight'*kron(M0Right,M1Right)*D10xRight/Jright*dydetaRight^2;	% contribution in x-dir.
        A2Right = D10yRight'*kron(M0Right,M1Right)*D10yRight/Jright*dxdxiRight^2; 	% contribution in y-dir.
        RHSright = kron(M0Right,M0Right)*fMeshRight*Jright;                        	% right hand side
        Aright = (A1Right+A2Right);

        % global system matrix
        A = blkdiag(Aleft,Aright);
        b = [RHSleft;RHSright];

        %-------------------------------------------------------------------------%
        % apply mortar method                                                     %
        %-------------------------------------------------------------------------%

        % lowest order
        if pLeft <= pRight
            pLow = pLeft;
        else
            pLow = pRight;
        end

        % count
        nInLeft = pLeft*(pLeft+1); 
        nInRight = pRight*(pRight+1);
        nMortar = pLow+1;
        nDOF = nInLeft+nInRight+nMortar;

        % boundary node 
        boundLocsLeft = (pLeft+1):(pLeft+1):(pLeft+1)^2;
        boundLocsRight = 1:(pRight+1):((pRight+1)^2-pRight);

        % inner nodes
        innerLocsLeft = 1:(pLeft+1)^2;
        innerLocsLeft(boundLocsLeft) = [];
        innerLocsRight = 1:(pRight+1)^2;
        innerLocsRight(boundLocsRight) = []; 

        % left element - Z-matrices
        ZiLeft = eye(nInLeft);
        WbLeft = Wleft; WmLeft = WbLeft*MimeticpolyVal(yMeshLeft,pLow,1)';
        ZbLeft = WbLeft\WmLeft;

        % right element - Z-matrices
        ZiRight = eye(nInRight);
        WbRight = Wright; WmRight = WbRight*MimeticpolyVal(yMeshRight,pLow,1)';
        ZbRight = WbRight\WmRight;

        % boundary nodes (of the inner nodes)
        bcLeft = unique( [1:pLeft ...                     	% lower side
                          1:pLeft:(1+pLeft^2) ...           % left side
                          (pLeft^2+1):(pLeft*(pLeft+1))] );	% upper side
        bcRight = unique( [1:pRight ...                             % lower side
                           pRight:pRight:(pRight*(pRight+1)) ...	% right side
                           pRight^2+1:(pRight*(pRight+1))] );   	% upper side

        % fill numbering matrices
        innerLeft = zeros(nInLeft,2);
        innerRight = zeros(nInRight,2);
        innerLeft(:,1) = innerLocsLeft;
        innerLeft(:,2) = 1:nInLeft;    
        innerRight(:,1) = innerLocsRight+nNodesLeft;
        innerRight(:,2) = (1:nInRight)+nInLeft;

        mortarC = (1:nMortar)+max(innerRight(:,2));       

        Zglobal = zeros(nNodesLeft+nNodesRight,nDOF);
        
        Zglobal(innerLeft(:,1),innerLeft(:,2)) = ZiLeft;
        Zglobal(innerRight(:,1),innerRight(:,2)) = ZiRight;
        Zglobal(boundLocsLeft,mortarC) = ZbLeft;
        Zglobal(boundLocsRight+nNodesLeft,mortarC) = ZbRight;            

        bcOut = [innerLeft(bcLeft,2); innerRight(bcRight,2); mortarC([1,end])'];

        %-------------------------------------------------------------------------%
        % matrix + RHS                                                            %
        %-------------------------------------------------------------------------%

        %%% ZAZ-matrix
        ZAZ = Zglobal'*A*Zglobal;
        Zb = Zglobal'*b;

        %-------------------------------------------------------------------------%
        % apply boundary conditions 	                                          %
        %-------------------------------------------------------------------------%        

        bcValues = zeros(length(bcOut),1);
        Zb = Zb-ZAZ(:,bcOut)*bcValues;
        ZAZ(bcOut,:) = [];                   
        ZAZ(:,bcOut) = [];                    
        Zb(bcOut) = [];   

        %-------------------------------------------------------------------------%
        % solve                                                                   %
        %-------------------------------------------------------------------------%

        ZAZ = sparse(ZAZ);
        Zb = sparse(Zb);
        x = ZAZ\Zb;
        % disp(x)

        %-------------------------------------------------------------------------%
        % reconstruct solution vector                                             %
        %-------------------------------------------------------------------------%

        nodesIn = 1:nDOF; nodesIn(bcOut) = [];
        sol = zeros(nDOF,1); 
        sol(nodesIn) = x;
        sol(bcOut) = bcValues;        

        totalSol = Zglobal*sol;
        solLeft = totalSol(1:nNodesLeft);
        solRight = totalSol((1:nNodesRight)+nNodesLeft);


        %-------------------------------------------------------------------------%
        % error determination                                                     %
        %-------------------------------------------------------------------------%    

        % left element error
        elementSolLeft = reshape(solLeft,pLeft+1,pLeft+1)';   
        errorLeft = normsTwoD(elementSolLeft,funcPar,0,mapMethodLeft,mapVarLeft,errorQuad);     

        % right element error
        elementSolRight = reshape(solRight,pRight+1,pRight+1)';   
        errorRight = normsTwoD(elementSolRight,funcPar,0,mapMethodRight,mapVarRight,errorQuad);  

        % global error
        errorGlobal = errorLeft+errorRight;

        %-------------------------------------------------------------------------%
        % print some output                                                       %
        %-------------------------------------------------------------------------%    

        fm = '%s %2.0f \t  %s %2.0f \t  %s %1.3e \t  %s %1.3e \t  %s %1.3e \n';
        fprintf(fm,'pL =',pLeft,'pR =',pRight,...
                'eL =',errorLeft,'eR =',errorRight,'e =',errorGlobal);
            
        %-------------------------------------------------------------------------%
        % plot results                                                            %
        %-------------------------------------------------------------------------%        
        %%

        plotResultsFlag = 1;
        if plotResultsFlag == 1

            %%% plotting the exact solution    
            figure(figureCount),subplot(1,3,1), hold on                                                           
            surf(xFine,yFine,U_Fine)
            shading flat, axis tight, grid on, axis equal,
            plotExactStr{1} = 'exact';
            plotExactStr{2} = strcat('fineGridN=',num2str(fineGridN));
            title(plotExactStr), xlabel('x'), ylabel('y'), zlabel('z'), view(3)

            %%% plot element corner nodes 
            subplot(1,3,2), hold on     
            % left element
            xEl = xElementCoordsLeft;
            yEl = yElementCoordsLeft;
            Z_ElementValues = funcTwoD(xEl,yEl,0,funcPar,1);
            [X_ElementCoords Y_ElementCoords] = meshgrid(xEl,yEl);
            plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2), view(3)       
            % right element
            xEl = xElementCoordsRight;
            yEl = yElementCoordsRight;
            Z_ElementValues = funcTwoD(xEl,yEl,0,funcPar,1);
            [X_ElementCoords Y_ElementCoords] = meshgrid(xEl,yEl);
            plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2), view(3)           

            %%% construct solution matrix
            % left element
            elementSol = reshape(solLeft,pLeft+1,pLeft+1)';
            interpSolLeft = hfLeft'*elementSol*hfLeft;            
            % right element
            elementSol = reshape(solRight,pRight+1,pRight+1)';
            interpSolRight = hfRight'*elementSol*hfRight;                

            %%% coloring
            myColor = zeros(length(totalSol),1);
            vmin = min(totalSol); vmax = max(totalSol); vrange = vmax-vmin';
            for k=1:length(totalSol)
                myColor(k) = (totalSol(k)-vmin)/vrange+1;
            end

            % left element
            [X_MeshZeroForms Y_MeshZeroForms] = ZeroFormCoordsTwoD(1,pLeft,xElementCoordsLeft,yElementCoordsLeft); 
            xScat = X_MeshZeroForms';
            yScat = Y_MeshZeroForms';                            
            scatter3(xScat,yScat,solLeft,25,myColor(1:nNodesLeft),'filled');
            % right element
            [X_MeshZeroForms Y_MeshZeroForms] = ZeroFormCoordsTwoD(1,pRight,xElementCoordsRight,yElementCoordsRight); 
            xScat = X_MeshZeroForms';
            yScat = Y_MeshZeroForms';                            
            scatter3(xScat,yScat,solRight,25,myColor((nNodesLeft+1):end),'filled');   
            axis tight, grid on, axis equal, title('nodal'), view(3), xlabel('x'), ylabel('y'), zlabel('z')

            %%% plot interpolated solution
            subplot(1,3,3), hold on
            [xFineEl yFineEl] = mappingMethod(mapMethodLeft,mapVarLeft,fineGridRef,fineGridRef); 
            surf(xFineEl,yFineEl,interpSolLeft); 
            [xFineEl yFineEl] = mappingMethod(mapMethodRight,mapVarRight,fineGridRef,fineGridRef); 
            surf(xFineEl,yFineEl,interpSolRight); 
            shading flat, axis tight, grid on, axis equal, view(3)                 

            figureCount = figureCount+1;                    
        end              
            
            
    end
end

    
