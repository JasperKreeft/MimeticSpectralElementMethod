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
pLeft = 2;  % left element order
pRight = 4; % right element order

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

funcPar = [1 1 1];	% [function, function parameters wx/wy]
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

% 0-forms
figure(figureCount)
% subplot(1,2,[1 2])
subplot(1,2,1), hold on
PlotGlobalNumberingZeroFormTwoD(hLeft,pLeft,xElementCoordsLeft,yElementCoordsLeft,1,1)
subplot(1,2,2), hold on
PlotGlobalNumberingZeroFormTwoD(hRight,pRight,xElementCoordsRight,yElementCoordsRight,1,1)
figureCount = figureCount+1;

% % 1-forms
% figure(figureCount)
% subplot(1,2,1), hold on
% PlotGlobalNumberingOneFormTwoD(hLeft,pLeft,xElementCoordsLeft,yElementCoordsLeft,[1 1],1)
% subplot(1,2,2), hold on
% PlotGlobalNumberingOneFormTwoD(hRight,pRight,xElementCoordsRight,yElementCoordsRight,[1 1],1)
% figureCount = figureCount+1;

% % 2-forms
% figure(figureCount)
% subplot(1,2,1), hold on
% PlotGlobalNumberingTwoFormTwoD(hLeft,pLeft,xElementCoordsLeft,yElementCoordsLeft,1,1)
% subplot(1,2,2), hold on
% PlotGlobalNumberingTwoFormTwoD(hRight,pRight,xElementCoordsRight,yElementCoordsRight,1,1)
% figureCount = figureCount+1;

%-------------------------------------------------------------------------%
% calculate exact solution (fine grid)                                    %
%-------------------------------------------------------------------------%

fineGridRef = linspace(-1,1,fineGridN);                                     % finegrid    
[xFine yFine] = mappingMethod(mapMethod,mapVar,fineGridRef,fineGridRef);	% interpolation mesh x-dir./y-dir.
U_Fine = funcTwoD(xFine,yFine,0,funcPar,1);                                 % exact solution (matrix) 
F_Fine = -funcTwoD(xFine,yFine,2,funcPar,1);

%-------------------------------------------------------------------------%
% allocate storage space                                                  %
%-------------------------------------------------------------------------%

...

%-------------------------------------------------------------------------%
% set up and solve systems                                                %
%-------------------------------------------------------------------------%

...

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
% mesh geometry                                                           %
%-------------------------------------------------------------------------%

...
    
%-------------------------------------------------------------------------%
% mesh topology                                                           % 
%-------------------------------------------------------------------------%    

%%% left element
D10xLeft = dZeroXiTwoD(pLeft);	% incidence matrix x-dir. 
D10yLeft = dZeroEtaTwoD(pLeft);	% incidence matrix y-dir.

ZeroFormsLeft = GlobalNumberingZeroFormTwoD(hLeft,pLeft);	% 0-forms local numbering
OneFormsLeft  = GlobalNumberingOneFormTwoD(hLeft,pLeft);	% 1-forms local numbering
TwoFormsLeft  = GlobalNumberingTwoFormTwoD(hLeft,pLeft);  	% 2-forms local numbering    

nNodesLeft = max(max(ZeroFormsLeft));	% # of 0-forms
nLinesLeft = max(max(OneFormsLeft));   	% # of 1-forms
nSurfacesLeft = max(max(TwoFormsLeft));	% # of 2-forms

%%% right element
D10xRight = dZeroXiTwoD(pRight);	% incidence matrix x-dir.
D10yRight = dZeroEtaTwoD(pRight);	% incidence matrix y-dir.

ZeroFormsRight = GlobalNumberingZeroFormTwoD(hRight,pRight);	% 0-forms local numbering
OneFormsRight  = GlobalNumberingOneFormTwoD(hRight,pRight);     % 1-forms local numbering
TwoFormsRight  = GlobalNumberingTwoFormTwoD(hRight,pRight);  	% 2-forms local numbering    

nNodesRight = max(max(ZeroFormsRight));     % # of 0-forms
nLinesRight = max(max(OneFormsRight));   	% # of 1-forms
nSurfacesRight = max(max(TwoFormsRight));	% # of 2-forms

%-------------------------------------------------------------------------%
% storage                                                                 %
%-------------------------------------------------------------------------%

boundaryNodesLnLeft = cell(3,1);    % local numbering boundary nodes left element 
boundaryNodesLnRight = cell(3,1);   % local numbering boundary nodes right element

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

%-------------------------------------------------------------------------%
% apply mortar method                                                     %
%-------------------------------------------------------------------------%

% lowest order
if pLeft <= pRight
    pLow = pLeft;
else
    pLow = pRight;
end

% boundary node 
boundLocsLeft = (pLeft+1):(pLeft+1):(pLeft+1)^2;
boundLocsRight = 1:(pRight+1):((pRight+1)^2-pRight);

% inner nodes
innerLocsLeft = 1:(pLeft+1)^2;
innerLocsLeft(boundLocsLeft) = [];
innerLocsRight = 1:(pRight+1)^2;
innerLocsRight(boundLocsRight) = []; 

% mortar nodes
if pLeft == pRight % elements are of the same order
    inCommonHigh = 1:pLeft+1;
    inCommonLow = inCommonHigh;
    mortarLocsLeft = boundLocsLeft;
    mortarLocsRight = boundLocsRight;
else % elements are not of the same order    
    if pRight>pLeft % higher order on the right          
        remainder = rem(pRight,pLeft);
        if remainder == 0 % pRight is integer multiple of pLeft
            inCommonHigh = 1:(pRight/pLeft):(pRight+1);
            inCommonLow = 1:(pLeft+1);
            mortarLocsLeft = boundLocsLeft;
            mortarLocsRight = boundLocsRight(inCommonHigh);
        else % pRight is NOT integer multiple of pLeft
            inCommonHigh = [1,pRight+1];
            inCommonLow = [1,pLeft+1];
            mortarLocsLeft = boundLocsLeft;
            mortarLocsRight = boundLocsRight(inCommonHigh);
        end    
    elseif pLeft>pRight % higher order on the left
        remainder = rem(pLeft,pRight);
        if remainder == 0 % pLeft is integer multiple of pRight
            inCommonHigh = 1:(pLeft/pRight):(pLeft+1);
            inCommonLow = 1:pRight+1;
            mortarLocsLeft = boundLocsLeft(inCommonHigh);
            mortarLocsRight = boundLocsRight;
        else % pLeft is NOT integer multiple of pRight
            inCommonHigh = [1,pLeft+1];
            inCommonLow = [1,pRight+1];
            mortarLocsLeft = boundLocsLeft(inCommonHigh);
            mortarLocsRight = boundLocsRight;
        end          
    end        
end
notInCommonHigh = 1:inCommonHigh(end);
notInCommonHigh(inCommonHigh) = [];
% notInCommonLow = 1:inCommonLow(end);%not used
% notInCommonLow(inCommonLow) = [];%not used

% left element - Z-matrices
ZiLeft = eye(nNodesLeft-(pLeft+1));
WbLeft = Wleft; WmLeft = WbLeft*MimeticpolyVal(yMeshLeft,pLow,1)';
ZbLeft = WbLeft\WmLeft;
% Zleft = blkdiag(ZiLeft,ZbLeft);

% right element - Z-matrices
ZiRight = eye(nNodesRight-(pRight+1));
WbRight = Wright; WmRight = WbRight*MimeticpolyVal(yMeshRight,pLow,1)';
ZbRight = WbRight\WmRight;
% Zright = blkdiag(ZiRight,ZbRight);

% system storage
A = zeros(nNodesLeft+nNodesRight-uint32(length(inCommonHigh)));	
b = zeros(nNodesLeft+nNodesRight-uint32(length(inCommonHigh)),1);     

% construct global Z and A
Zglobal = zeros((nNodesLeft+nNodesRight),(nNodesLeft+nNodesRight));
% identities of inner nodes (left and right)
Zglobal(innerLocsLeft,innerLocsLeft) = ZiLeft;
Zglobal(uint32(innerLocsRight)+nNodesLeft,uint32(innerLocsRight)+nNodesLeft) = ZiRight;
% adding boundary nodes (left and right)
if pLeft<=pRight
    
    %%% Z-matrix
    % adding boundary nodes of left element (identity)
    Zglobal(mortarLocsLeft,mortarLocsLeft) = ZbLeft;
    % adding boundary nodes of right element - mortar nodes (the ones which it has in common, always at least the corners)
    Zglobal(boundLocsLeft(inCommonLow),boundLocsLeft(inCommonLow))= ...
        Zglobal(boundLocsLeft(inCommonLow),boundLocsLeft(inCommonLow))+eye(length(inCommonLow));
    % adding boundary nodes of right element - remaining nodes
    Zglobal(uint32(boundLocsRight(notInCommonHigh))+nNodesLeft,mortarLocsLeft)= ZbRight(notInCommonHigh,:);
    % determine which columns and rows to remove
    removeColumns = nNodesLeft+uint32(boundLocsRight);%remove all boundary nodes of the higher order element
    removeRows = nNodesLeft+uint32(mortarLocsRight);%remove the duplicate mortar nodes
    % final Z-matrix    
    Zglobal(:,removeColumns) = [];
    Zglobal(removeRows,:) = [];
       
    %%% global numbering
    globalLeft = 1:nNodesLeft;    
    globalRight = zeros(1,nNodesRight);
    globalRight(mortarLocsRight) = boundLocsLeft(inCommonLow);
    uniqueRight = 1:nNodesRight;
    uniqueRight(mortarLocsRight) = [];
    globalRight(uniqueRight) = (1:(nNodesRight-length(inCommonHigh)))+nNodesLeft;
    
    % explanation: 
    % mortarLocsLeft = locations of the mortar nodes of the left element
    % mortarLocsRight = locations of the mortar nodes of the right element
    % boundLocsLeft = locations of the boundary nodes of the left element
    % boundLocsRight = locations of the boundary nodes of the right element
    % inCommonLow = mortar nodes which the left element has in common with the right element (numbers along boundary nodes)
    % inCommonHigh = mortar nodes which the right element has in common with the left element (numbers along boundary nodes)
    % notInCommonLow = NOT USED/DEFINED
    % notInCommonHigh = boundary nodes of the left element which are not mortar nodes
    % use: apply inCommon(Low/High) and notInCommonHigh to boundLocs(Left/Right)
    % note: boundLocsLeft(inCommonLow) ~= mortarLocsLeft
    %       boundLocsRight(inCommonHigh) = mortarLocsRight
    
%     %%% ordering
%     newOrderLeft = 1:nNodesLeft;    
%     newOrderRight = (1:(nNodesRight-length(boundLocsRight)))+nNodesLeft;
        
    %%% boundary nodes 
    % left element
    boundaryNodesLnLeft{1} = 1:(pLeft+1);% lower boundary
    boundaryNodesLnLeft{2} = 1:(pLeft+1):(1+pLeft*(pLeft+1));% left boundary
    boundaryNodesLnLeft{3} = boundaryNodesLnLeft{1}+(pLeft*(pLeft+1));% upper boundary
    bcLeft = unique([boundaryNodesLnLeft{1} boundaryNodesLnLeft{2}  boundaryNodesLnLeft{3}]);
    lastVal = max(bcLeft);
    % right element 
%     boundaryNodesLnRight{1} = globalRight(2:(pRight+1));% lower boundary
%     boundaryNodesLnRight{2} = globalRight( ((pRight+1):(pRight+1):(pRight+1)^2) );% right boundary
%     boundaryNodesLnRight{3} = globalRight( (pRight*(pRight+1)+2):(pRight+1)^2 );% upper boundary    
    boundaryNodesLnRight{1} = lastVal+(1:pRight);% lower boundary
    boundaryNodesLnRight{2} = lastVal+(pRight:pRight:(pRight+1)*pRight);% right boundary
    boundaryNodesLnRight{3} = lastVal+((pRight^2+1):(pRight*(pRight+1)));% upper boundary
    bcRight = unique([boundaryNodesLnRight{1} boundaryNodesLnRight{2}  boundaryNodesLnRight{3}]);
    
    
    
% bcLeft = unique( [1:(pLeft+1) ...                             % lower side
%                   1:(pLeft+1):(1+pLeft*(pLeft+1)) ...         % left side
%                   (1+pLeft*(pLeft+1)):((pLeft+1)^2)] );       % upper side
% bcRight = unique( [1:(pRight+1) ...                           % lower side
%                    (pRight+1):(pRight+1):((pRight+1)^2) ...	% right side
%                    (1+pRight*(pRight+1)):((pRight+1)^2)] );   % upper side
    
else % pLeft>pRight
    
    %%% Z-matrix
    % adding boundary nodes of right element (identity)
    Zglobal(uint32(mortarLocsRight)+nNodesLeft,uint32(mortarLocsRight)+nNodesLeft) = ZbRight;    
    % adding boundary nodes of left element - mortar nodes (the ones which it has in common, always at least the corners)
    Zglobal(uint32(mortarLocsRight)+nNodesLeft,uint32(mortarLocsRight)+nNodesLeft)= ...
        Zglobal(uint32(mortarLocsRight)+nNodesLeft,uint32(mortarLocsRight)+nNodesLeft)+eye(length(mortarLocsRight));      
    % adding boundary nodes of left element - remaining nodes
    Zglobal(boundLocsLeft(notInCommonHigh),uint32(mortarLocsRight)+nNodesLeft)= ZbLeft(notInCommonHigh,:);
    % remove rows and columns
    removeColumns = boundLocsLeft;%remove all boundary nodes of the higher order element
    removeRows = mortarLocsLeft;%remove the duplicate mortar nodes  
    % final Z-matrix    
    Zglobal(:,removeColumns) = [];
    Zglobal(removeRows,:) = [];
    
    %%% global numbering
    globalRight = (1:nNodesRight)+(nNodesLeft-length(inCommonHigh));
    globalLeft = zeros(1,nNodesLeft);    
    globalLeft(mortarLocsLeft) = globalRight(inCommonLow);    
    uniqueLeft = 1:nNodesLeft;
    uniqueLeft(mortarLocsLeft) = [];   
    globalLeft(uniqueLeft) = 1:(nNodesLeft-length(inCommonHigh));    
    
    %%% ordering
    newOrderLeft = 1:(nNodesLeft-length(boundLocsLeft));    
    newOrderRight = (1:nNodesRight)+newOrderLeft(end);
   
    %%% boundary nodes
    % left element
    boundaryNodesLnLeft{1} = 1:pLeft;% lower boundary
    boundaryNodesLnLeft{2} = 1:pLeft:(pLeft*(pLeft+1));% left boundary
    boundaryNodesLnLeft{3} = boundaryNodesLnLeft{1}+(pLeft^2);% upper boundary
    bcLeft = unique([boundaryNodesLnLeft{1} boundaryNodesLnLeft{2}  boundaryNodesLnLeft{3}]);
    lastVal = max(bcLeft);
    % right element
    boundaryNodesLnRight{1} = lastVal+(1:(pRight+1));% lower boundary
    boundaryNodesLnRight{2} = lastVal+((pRight+1):(pRight+1):(pRight+1)^2);% right boundary    
    boundaryNodesLnRight{3} = lastVal+((pRight*(pRight+1)+1):(pRight+1)^2);% upper boundary
    bcRight = unique([boundaryNodesLnRight{1} boundaryNodesLnRight{2}  boundaryNodesLnRight{3}]);

end

%-------------------------------------------------------------------------%
% global numbering                                                        %
%-------------------------------------------------------------------------%

...

%-------------------------------------------------------------------------%
% matrix + RHS                                                            %
%-------------------------------------------------------------------------%

%%% A-matrix (and RHS)    
A(globalLeft,globalLeft) =  A(globalLeft,globalLeft)+Aleft;
A(globalRight,globalRight) =  A(globalRight,globalRight)+Aright;
b(globalLeft,1) = b(globalLeft,1)+RHSleft;
b(globalRight,1) = b(globalRight,1)+RHSright;

%%% ZAZ-matrix
ZAZ = Zglobal'*A*Zglobal;
Zb = Zglobal'*b;

%%% total number of nodes to be calculated
nNodes = length(Zb);

%-------------------------------------------------------------------------%
% numbering of boundary nodes/elements                                    %
%-------------------------------------------------------------------------%

...
 
%---------------------------------------------------------%
% apply boundary conditions 	                          %
%---------------------------------------------------------%        

bcValues = zeros(length(bcLeft)+length(bcRight),1);
totalBC = [bcLeft bcRight];

Zb = Zb-ZAZ(:,totalBC)*bcValues;
ZAZ(totalBC,:) = [];                   
ZAZ(:,totalBC) = [];                    
Zb(totalBC) = [];   

%---------------------------------------------------------%
% solve                                                   %
%---------------------------------------------------------%

ZAZ = sparse(ZAZ);
Zb = sparse(Zb);
x = ZAZ\Zb;
% disp(x)

%---------------------------------------------------------%
% reconstruct solution vector                             %
%---------------------------------------------------------%    

nodesIn = 1:nNodes; nodesIn(totalBC) = [];
sol = zeros(nNodes,1); 
sol(nodesIn) = x;
sol(totalBC) = bcValues;        

if pLeft<=pRight
    sol2 = Zglobal*sol;
    
    sol2left = sol2(1:nNodesLeft);
    sol2right = sol2((nNodesLeft+1):end);% does not include mortar nodes (the ones inCommonHigh)
    
    solLeft = sol2left;     
    solRight = zeros(nNodesRight,1);
    
    test = 1:nNodesRight;
    test(mortarLocsRight) = [];
    mortarVals = solLeft(boundLocsLeft(inCommonLow));
    solRight(test) = sol2right;
    solRight(mortarLocsRight) = mortarVals;
    solRight(boundLocsRight(notInCommonHigh))
%     solRight = zeros(nNodesRight,1);
%     solRight(innerLocsRight) = sol2((nNodesLeft+1):end);
%     solRight(mortarLocsRight) = solLeft(mortarLocsLeft);
%     solLeft = sol(1:nNodesLeft);
%     solRight = zeros(nNodesRight,1);
%     solRight(innerLocsRight) = sol((nNodesLeft+1):end);
%     % add mortar solution
%     mortarSol = ZbRight*solLeft(mortarLocsLeft);
%     solRight(boundLocsRight) = mortarSol;    
else % pLeft>pRight 
    solRight = sol((nNodesLeft-length(boundLocsLeft)+1):end);
    solLeft = zeros(nNodesLeft,1);
    solLeft(innerLocsLeft) = sol(1:(nNodesLeft-length(boundLocsLeft)));
    % add mortar solution
    mortarSol = ZbLeft*solRight(mortarLocsRight);
    solLeft(boundLocsLeft) = mortarSol;    
end
solGlobal = [solLeft;solRight];

%---------------------------------------------------------%
% error determination                                     %
%---------------------------------------------------------%    

...

%---------------------------------------------------------%
% print output                                            %
%---------------------------------------------------------%        

...             

%---------------------------------------------------------%
% plot results                                            %
%---------------------------------------------------------%        
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
    elementSol = solLeft;
    elementSol = reshape(elementSol,pLeft+1,pLeft+1)';
    interpSolLeft = hfLeft'*elementSol*hfLeft;            
    % right element
    elementSol = solRight;
    elementSol = reshape(elementSol,pRight+1,pRight+1)';
    interpSolRight = hfRight'*elementSol*hfRight;                
        
    %%% coloring
    myColor = zeros(length(solGlobal),1);
    vmin = min(solGlobal); vmax = max(solGlobal); vrange = vmax-vmin';
    for k=1:length(solGlobal)
        myColor(k) = (solGlobal(k)-vmin)/vrange+1;
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

%         plotNodalStr{1} = 'nodal';
%         plotNodalStr{2} = strcat('h_x=',num2str(hx),', h_y=',num2str(hy));                            
%         plotNodalStr{3} = strcat('p_x=',num2str(px),', p_y=',num2str(py));
%         plotNodalStr{4} = strcat('BC=[',num2str(BCs),']');
%         title(plotNodalStr), 

    %%% fine coordinates

    %%% plot interpolated solution
    subplot(1,3,3), hold on
    [xFineEl yFineEl] = mappingMethod(mapMethodLeft,mapVarLeft,fineGridRef,fineGridRef); 
    surf(xFineEl,yFineEl,interpSolLeft); 
    [xFineEl yFineEl] = mappingMethod(mapMethodRight,mapVarRight,fineGridRef,fineGridRef); 
    surf(xFineEl,yFineEl,interpSolRight); 
%     plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2)
    shading flat, axis tight, grid on, axis equal, view(3)
    
%     plotInterpStr{1} = 'interpolated';
%     plotInterpStr{2} = strcat('h_x=',num2str(hx),', h_y=',num2str(hy));                            
%     plotInterpStr{3} = strcat('p_x=',num2str(px),', p_y=',num2str(py));
%     plotInterpStr{4} = strcat('BC=[',num2str(BCs),']');            
%     title(plotInterpStr), xlabel('x'), ylabel('y'), zlabel('z')
% 
%     %%% plot element boundaries
%     X_FineV = repmat(xElementCoords',fineGridN,1);
%     Y_FineV = repmat(yFine,1,hx+1);
%     Z_FineV = funcTwoD(X_FineV,Y_FineV,0,funcPar,3);
%     plot3(X_FineV,Y_FineV,Z_FineV,'k--','linewidth',3)
% 
%     Y_FineH = repmat(yElementCoords',fineGridN,1);
%     X_FineH = repmat(xFine,1,hy+1);
%     Z_FineH = funcTwoD(X_FineH,Y_FineH,0,funcPar,3);
%     plot3(X_FineH,Y_FineH,Z_FineH,'k--','linewidth',3)                              
                
    figureCount = figureCount+1;                    
end      
