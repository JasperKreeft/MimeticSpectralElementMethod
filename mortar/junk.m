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
pRight = 2; % right element order

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
subplot(1,2,[1 2])
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
[hnleft enleft] = MimeticpolyVal(refMeshNodesLeft,pLeft,1);	% coarse mesh interpolants (x-dir.) 
[hfleft efleft] = MimeticpolyVal(fineGridRef,pLeft,1);      % fine mesh interpolants (x-dir.)

%%% right element
[refMeshNodesRight refMeshWeightsRight] = GLLnodes(pRight);     % mesh reference nodes and weights (x-dir.)
[hnright enright] = MimeticpolyVal(refMeshNodesRight,pRight,1);	% coarse mesh interpolants (x-dir.) 
[hfright efright] = MimeticpolyVal(fineGridRef,pRight,1);       % fine mesh interpolants (x-dir.)

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

% check for duplicate nodes/lines
if pLeft == pRight % elements are of the same order   
    duplicateNodes = pLeft+1;
    % duplicateNodes = pRight+1;
    duplicateLines = pLeft;
    % duplicateLines = pRight;
else % elements are not of the same order    
    if pRight>pLeft % higher order on the right          
        remainder = rem(pRight,pLeft);
        if remainder == 0 % pRight is integer multiple of pLeft
            duplicateNodes = pLeft+1;
            duplicateLines = pLeft;
        else % pRight is NOT integer multiple of pLeft
            % all unique, except for 2 corner end nodes
            duplicateNodes = 2;
            duplicateLines = -(pLeft-1); % it actually creates extra line segments, hence it's <0
        end    
    elseif pLeft>pRight % higher order on the left
        remainder = rem(pLeft,pRight);
        if remainder == 0 % pLeft is integer multiple of pRight
            duplicateNodes = pRight+1;
            duplicateLines = pRight;
        else % pLeft is NOT integer multiple of pRight
            % all unique, except for 2 corner end nodes
            duplicateNodes = 2;
            duplicateLines = -(pRight-1); % it actually creates extra line segments, hence it's <0
        end          
    end        
end

% global values
nNodes = nNodesLeft+nNodesRight-duplicateNodes;
nLines = nLinesLeft+nLinesRight-duplicateLines;
nSurfaces = nSurfacesLeft+nSurfacesRight;

% % check for duplicate nodes/lines
% if pLeft == pRight % elements are of the same order   
%     duplicateNodes = pLeft+1;
%     % duplicateNodes = pRight+1;
%     duplicateLines = pLeft;
%     % duplicateLines = pRight;
% else % elements are not of the same order    
%     if pRight>pLeft % higher order on the right          
%         remainder = rem(pRight,pLeft);
%         if remainder == 0 % pRight is integer multiple of pLeft
%             duplicateNodes = pLeft+1;
%             duplicateLines = pLeft;            
%         else % pRight is NOT integer multiple of pLeft
%             % all unique, except for 2 corner end nodes
%             duplicateNodes = 2;
%             duplicateLines = -(pLeft-1); % it actually creates extra line segments, hence it's <0                  
%         end    
%     elseif pLeft>pRight % higher order on the left
%         remainder = rem(pLeft,pRight);
%         if remainder == 0 % pLeft is integer multiple of pRight
%             duplicateNodes = pRight+1;
%             duplicateLines = pRight;                
%         else % pLeft is NOT integer multiple of pRight
%             % all unique, except for 2 corner end nodes
%             duplicateNodes = 2;
%             duplicateLines = -(pRight-1); % it actually creates extra line segments, hence it's <0                  
%         end          
%     end        
% end

% % global values
% nNodes = nNodesLeft+nNodesRight-duplicateNodes;
% nLines = nLinesLeft+nLinesRight-duplicateLines;
% nSurfaces = nSurfacesLeft+nSurfacesRight;

% % boundary nodes (local numbering)
% leftElVal = (pLeft+1):(pLeft+1):(pLeft+1)^2;
% leftElValMortar = leftElVal;
% rightElVal = 1:(pRight+1):((pRight+1)^2-pRight); 
% RightElValMortar = rightElVal;

%-------------------------------------------------------------------------%
% storage                                                                 %
%-------------------------------------------------------------------------%

...

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
dx_dxiLeft = Dleft(1,1);                           
dy_detaLeft = Dleft(2,2);
% matrix construction
Wleft = spdiags(refMeshWeightsLeft,0,pLeft+1,pLeft+1);                  % integration weights x-dir.
M0left = hnleft*Wleft*hnleft';                                          % lagrange integration x-dir.
M1left = enleft*Wleft*enleft';                                          % edge integration x-dir.
A1left = D10xLeft'*kron(M0left,M1left)*D10xLeft/Jleft*dy_detaLeft^2;	% contribution in x-dir.
A2left = D10yLeft'*kron(M0left,M1left)*D10yLeft/Jleft*dx_dxiLeft^2;     % contribution in y-dir.
RHSleft = kron(M0left,M0left)*fMeshLeft*Jleft;                          % right hand side

%%% right element
% mesh
[xMeshRight yMeshRight] = mappingMethod(mapMethodRight,mapVarRight,refMeshNodesRight,refMeshNodesRight);
% RHS
fMeshRight = -funcTwoD(xMeshRight,yMeshRight,2,funcPar,2);
% jacobian                           
[Jright Dright] = jacobianTwoD(mapMethodRight,mapVarRight);	% NOTE: J is just a constant for 'scale'
dx_dxiRight = Dright(1,1);                           
dy_detaRight = Dright(2,2);  
% matrix construction
Wright = spdiags(refMeshWeightsRight,0,pRight+1,pRight+1);                  % integration weights x-dir.
M0right = hnright*Wright*hnright';                                          % lagrange integration x-dir.
M1right = enright*Wright*enright';                                          % edge integration x-dir.
A1right = D10xRight'*kron(M0right,M1right)*D10xRight/Jright*dy_detaRight^2;	% contribution in x-dir.
A2right = D10yRight'*kron(M0right,M1right)*D10yRight/Jright*dx_dxiRight^2; 	% contribution in y-dir.
RHSright = kron(M0right,M0right)*fMeshRight*Jright;                        	% right hand side

% % adjust A to Z
% leftOrdering = [ZeroFormsLeft(innerLocsLeft) ZeroFormsLeft(boundLocsLeft)];
% Aleft = Aleft(leftOrdering,leftOrdering);
% rightOrdering = [ZeroFormsRight(innerLocsRight) ZeroFormsRight(boundLocsRight)];
% Aleft = Aleft(leftOrdering,leftOrdering);

% % construct ZAZ
% ZAZleft = Zleft'*Aleft*Zleft;
% bLeft = Zleft'*RHSleft;    
% ZAZright = Zright'*Aright*Zright;         
% bRight = Zright'*RHSright;    

% check
% innerNodesLeft = length(innerLocsLeft);
% innerNodesRight = length(innerLocsRight);
% innerNodes = innerNodesLeft+innerNodesRight;
% mortarNodes = pLow+1;
% nNodes = innerNodes+mortarNodes;

%-------------------------------------------------------------------------%
% global numbering                                                        %
%-------------------------------------------------------------------------%

% innerValsLeft = ZeroFormsLeft(innerLocsLeft);
% innerValsRight = ZeroFormsRight(innerLocsRight);
% valsLeft = 1:innerNodesLeft;
% valsRight = (1:innerNodesRight)+valsLeft(end);
% mortars = (1:mortarNodes)+valsRight(end);

%-------------------------------------------------------------------------%
% matrix + RHS                                                            %
%-------------------------------------------------------------------------%

% % storage
% A = zeros(nNodes,nNodes);
% b = zeros(nNodes,1);
% 
% % innerleft
% A(1:innerNodesLeft,1:size(ZAZleft)) = A(1:innerNodesLeft,1:size(ZAZleft))+...
%     ZAZleft(1:innerNodesLeft,:);
% b(1:innerNodesLeft) = b(1:innerNodesLeft)+bLeft(1:innerNodesLeft);
% % innerright
% A((innerNodesLeft+1):innerNodes,(innerNodesLeft+1):innerNodes) = A((innerNodesLeft+1):innerNodes,(innerNodesLeft+1):innerNodes)+...
%     ZAZleft(1:innerNodesRight,1:innerNodesRight);
% b((innerNodesLeft+1):innerNodes) = b((innerNodesLeft+1):innerNodes)+bRight(1:innerNodesRight);
% % mortarleft
% 
% % mortarright
% 
% 
% % build global matrix
% A(ZeroLeft,ZeroLeft) = A(ZeroLeft,ZeroLeft)+Aleft;
% A(ZeroRight,ZeroRight) = A(ZeroRight,ZeroRight)+Aright;
% b(ZeroLeft) = b(ZeroLeft)+bLeft;
% b(ZeroRight) = b(ZeroRight)+bRight;

%-------------------------------------------------------------------------%
% global numbering                                                        %
%-------------------------------------------------------------------------%

...

%-------------------------------------------------------------------------%
% apply mortar method                                                     %
%-------------------------------------------------------------------------%

% % relate mortars to local values
% if pLeft == pRight % elements are of the same order
%     inCommon = 1:pLeft+1;
%     mortarLocsLeft = boundLocsLeft;
%     mortarLocsRight = boundLocsRight;
% else % elements are not of the same order    
%     if pRight>pLeft % higher order on the right          
%         remainder = rem(pRight,pLeft);
%         if remainder == 0 % pRight is integer multiple of pLeft
%             inCommon = 1:(pRight/pLeft):(pRight+1);
%             mortarLocsLeft = boundLocsLeft;
%             mortarLocsRight = boundLocsRight(inCommon);
%         else % pRight is NOT integer multiple of pLeft
%             inCommon = [1,pRight+1];
%             mortarLocsLeft = boundLocsLeft;
%             mortarLocsRight = boundLocsRight([1,pRight+1]);
%         end    
%     elseif pLeft>pRight % higher order on the left
%         remainder = rem(pLeft,pRight);
%         if remainder == 0 % pLeft is integer multiple of pRight
%             inCommon = 1:(pLeft/pRight):(pLeft+1);
%             mortarLocsLeft = boundLocsLeft(inCommon);
%             mortarLocsRight = boundLocsRight;
%         else % pLeft is NOT integer multiple of pRight
%             inCommon = [1,pLeft+1];
%             mortarLocsLeft = boundLocsLeft([1,pLeft+1]);
%             mortarLocsRight = boundLocsRight;
%         end          
%     end        
% end


% Zleft = zeros(size(blkdiag(ZiLeft,ZbLeft)));
% Zleft(internalLocsLeft,internalLocsLeft) = ZiLeft;
% Zleft(boundLocsLeft,mortarLocsLeft) = ZbLeft;
% logicalLeft = any(Zleft);
% Zleft = Zleft(:,logicalLeft);

% Zright = zeros(size(blkdiag(ZiRight,ZbRight)));
% Zright(internalLocsRight,internalLocsRight) = ZiRight;
% Zright(boundLocsRight,mortarLocsRight) = ZbRight;
% logicalRight = any(Zright);
% Zright = Zright(:,logicalRight);

%%% left element
% boundary numbers (local)
leftElVal = (pLeft+1):(pLeft+1):(pLeft+1)^2; % local numbering
Zi = eye(2);
Wb = Wleft; Wm = Wb*MimeticpolyVal(yMeshLeft,pLeft,1);
Zb = Wb\Wm;

%%% right element
% boundary numbers (local)
rightElVal = 1:(pRight+1):((pRight+1)^2-pRight); % local numbering
innerNodes = 1:(pRight+1)^2;
innerNodes(rightElVal) = [];
W = Wright;
hn = MimeticpolyVal(yMeshRight,pLeft,1);	% coarse mesh interpolants (x-dir.) 
Wwiggle = Wright*hn';
Zright = W\Wwiggle;
Zid = eye(6);
ZZright = blkdiag(Zid,Zright);
% ZZrightNew = zeros(9,8);
% ZZrightNew(innerNodes,:) = ZZright(1:length(innerNodes),:);
% ZZrightNew(rightElVal,:) = ZZright(length(innerNodes)+1:end,:);
% ZZrightNew(innerNodes,:) = ZZright(1:length(innerNodes),:);
% ZZrightNew(rightElVal,:) = ZZright(length(innerNodes)+1:end,:);
ZZrightNew = ...
[1 0 0 0 0 0 0 0;
 0 1 0 0 0 0 0 0;
 0 0 1 0 0 0 0 0;
 .5 0 0 0 0 .5 0 0;
 0 0 0 1 0 0 0 0;
 0 0 0 0 1 0 0 0;
 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 1 0;
 0 0 0 0 0 0 0 1]; 
%-------------------------------------------------------------------------%
% store in matrix                                                         %
%-------------------------------------------------------------------------%

% left element
Aleft = ZZleft'*(A1left+A2left)*ZZleft;
bLeft = ZZleft'*RHSleft;    
% Aleft(ZeroFormsLeft,ZeroFormsLeft) = ZZleft'*(A1left+A2left)*ZZleft;
% bLeft(ZeroFormsLeft,1) = ZZleft'*RHSleft;    

% right element
ZeroFormsRightTemp = ZeroFormsRight;
ZeroFormsRightTemp(rightElVal(2:end-1)) = [];
% correct subsequent index numbers for removing as done before
for i=1:length(ZeroFormsRightTemp)
    if ZeroFormsRightTemp(i) > ZeroFormsRight(rightElVal(2:end-1))
        ZeroFormsRightTemp(i) = ZeroFormsRightTemp(i)-1;
    end
end
rightElValTemp = rightElVal;
for i=1:length(rightElValTemp)
    if rightElValTemp(i) > rightElVal(2:end-1)
        rightElValTemp(i) = rightElValTemp(i)-1;
    end
end
% Aright = zeros(8,8);
Aright = ZZrightNew'*(A1right+A2right)*ZZrightNew;         
bRight = ZZrightNew'*RHSright;    
% Aright(ZeroFormsRightTemp,ZeroFormsRightTemp) = ZZrightNew'*(A1right+A2right)*ZZrightNew;         
% bRight(ZeroFormsRightTemp,1) = ZZrightNew'*RHSright;      

% now this must really be checked!!!
% global numbering * note that mid-point is left out (so 10 is removed, 11 becomes new 10)
ZeroLeft = [1 9 2 10];
ZeroRight = [9 3 4 5 6 10 7 8];

% build global system matrix
A = zeros(10,10);
b = zeros(10,1);
A(ZeroLeft,ZeroLeft) = A(ZeroLeft,ZeroLeft)+Aleft;
A(ZeroRight,ZeroRight) = A(ZeroRight,ZeroRight)+Aright;
b(ZeroLeft) = b(ZeroLeft)+bLeft;
b(ZeroRight) = b(ZeroRight)+bRight;

%-------------------------------------------------------------------------%
% numbering of boundary nodes/elements                                    %
%-------------------------------------------------------------------------%

% %%% local numbering 
% boundaryNodesLn{1} = 1:(px+1);                          % lower boundary
% boundaryNodesLn{2} = boundaryNodesLn{1}+(py*(px+1));	% upper boundary
% boundaryNodesLn{3} = 1:(px+1):(1+py*(px+1));            % left boundary
% boundaryNodesLn{4} = boundaryNodesLn{3}+px;             % right boundary
% 
% %%% boundary elements
% boundaryElements{1} = 1:hx;                             % lower boundary          
% boundaryElements{2} = boundaryElements{1}+(hy-1)*hx;    % upper boundary
% boundaryElements{3} = 1:hx:(1+hx*(hy-1));               % left boundary
% boundaryElements{4} = boundaryElements{3}+(hx-1);      	% right boundary          
% 
% %%% global numbering 
% boundaryNodesGn{1} = ZeroForms(boundaryElements{1},boundaryNodesLn{1})'; 	% lower boundary 
% boundaryNodesGn{2} = ZeroForms(boundaryElements{2},boundaryNodesLn{2})';	% upper boundary 
% boundaryNodesGn{3} = ZeroForms(boundaryElements{3},boundaryNodesLn{3})';	% left boundary
% boundaryNodesGn{4} = ZeroForms(boundaryElements{4},boundaryNodesLn{4})';	% right boundary              

%---------------------------------------------------------%
% apply boundary conditions 	                          %
%---------------------------------------------------------%        

toRemove = [1 9 2 10 3 4 6 7 8];
bcValues = zeros(length(toRemove),1);
b = b-A(:,toRemove)*bcValues;
A(toRemove,:) = [];                   
A(:,toRemove) = [];                    
b(toRemove) = [];   

%         %%% update system rhs
%         b = b-A(:,dirichletNodes)*dirichletValues(firstVal); 
% 
%         %%% update system matrix
%         A(dirichletNodes,:) = [];                   
%         A(:,dirichletNodes) = [];                    
%         b(dirichletNodes) = [];   

% %%% dirichlet
% for nD = dirichletSides           
% boundaryValues{nD} = ...
% funcTwoD(X_MeshZeroForms(boundaryElements{nD},boundaryNodesLn{nD}),...
% Y_MeshZeroForms(boundaryElements{nD},boundaryNodesLn{nD}),0,funcPar,3);
% end
% %%% neumann 
% for nN = neumannSides
% ...
% end
% %%% call bc function
% [A b nodesOut firstVal dirichletValues] = ...
% bcZeroFormTwoD(BCs,A,b,boundaryNodesGn,boundaryValues);

%---------------------------------------------------------%
% solve                                                   %
%---------------------------------------------------------%

% error(elementCountX,elementCountY,orderCountX,orderCountY,2) = cond(A);
A = sparse(A);
b = sparse(b);
x = A\b;

%---------------------------------------------------------%
% reconstruct solution vector                             %
%---------------------------------------------------------%    

% if dirichletFlag
% nodesIn = 1:nNodes; nodesIn(nodesOut) = [];
% sol = zeros(nNodes,1); 
% sol(nodesIn) = x;
% sol(nodesOut) = dirichletValues(firstVal);        
% else
% sol = x;
% end

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
    
    subplot(1,3,2), hold on                                              
    surf(xFine,yFine,F_Fine)
    shading flat, axis tight, grid on, axis equal,
    plotExactStr{1} = 'exact - derivative';
    plotExactStr{2} = strcat('fineGridN=',num2str(fineGridN));
    title(plotExactStr), xlabel('x'), ylabel('y'), zlabel('z'), view(3)

%     %%% loop over elements, interpolate and plot results
%     for elementX=1:hx
%         for elementY=1:hy
% 
%         %%% global element numbering
%         element = elementX+hx*(elementY-1);
% 
%         %%% plot element corner nodes
%         subplot(1,3,2), hold on        
%         xEl = [xElementCoords(elementX) xElementCoords(elementX+1)];
%         yEl = [yElementCoords(elementY) yElementCoords(elementY+1)];
%         Z_ElementValues = funcTwoD(xEl,yEl,0,funcPar,1);
%         [X_ElementCoords Y_ElementCoords] = meshgrid(xEl,yEl);
%         plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2), view(3)                                                         
% 
%         %%% construct solution matrix
%         elementSol = sol(ZeroForms(element,:));
%         elementSol = reshape(elementSol,px+1,py+1)';
%         interpSol = hfy'*elementSol*hfx;            
% 
%         %%% coloring
%         myColor = zeros(length(sol),1);
%         vmin = min(sol); vmax = max(sol); vrange = vmax-vmin';
%         for k=1:length(sol)
%             myColor(k) = (sol(k)-vmin)/vrange+1;
%         end
% 
%         %%% nodal coordinates   
%         xScat = X_MeshZeroForms(element,:)';
%         yScat = Y_MeshZeroForms(element,:)';                            
% 
%         %%% plot solution nodes
%         myscat = scatter3(xScat,yScat,sol(ZeroForms(element,:)),25,myColor(ZeroForms(element,:)),'filled');
%         plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2)
%         axis tight, grid on, axis equal, title('nodal'), view(3)
%         plotNodalStr{1} = 'nodal';
%         plotNodalStr{2} = strcat('h_x=',num2str(hx),', h_y=',num2str(hy));                            
%         plotNodalStr{3} = strcat('p_x=',num2str(px),', p_y=',num2str(py));
%         plotNodalStr{4} = strcat('BC=[',num2str(BCs),']');
%         title(plotNodalStr), xlabel('x'), ylabel('y'), zlabel('z')
% 
%         %%% fine coordinates
%         mapVarEl = [xEl yEl];
%         [xFineEl yFineEl] = mappingMethod(mapMethod,mapVarEl,fineGridRef,fineGridRef); 
% 
%         %%% plot interpolated solution
%         subplot(1,3,3), hold on
%         surf(xFineEl,yFineEl,interpSol); 
%         plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2)
%         shading flat, axis tight, grid on, axis equal, view(3)
%         plotInterpStr{1} = 'interpolated';
%         plotInterpStr{2} = strcat('h_x=',num2str(hx),', h_y=',num2str(hy));                            
%         plotInterpStr{3} = strcat('p_x=',num2str(px),', p_y=',num2str(py));
%         plotInterpStr{4} = strcat('BC=[',num2str(BCs),']');            
%         title(plotInterpStr), xlabel('x'), ylabel('y'), zlabel('z')
% 
%         %%% plot element boundaries
%         X_FineV = repmat(xElementCoords',fineGridN,1);
%         Y_FineV = repmat(yFine,1,hx+1);
%         Z_FineV = funcTwoD(X_FineV,Y_FineV,0,funcPar,3);
%         plot3(X_FineV,Y_FineV,Z_FineV,'k--','linewidth',3)
% 
%         Y_FineH = repmat(yElementCoords',fineGridN,1);
%         X_FineH = repmat(xFine,1,hy+1);
%         Z_FineH = funcTwoD(X_FineH,Y_FineH,0,funcPar,3);
%         plot3(X_FineH,Y_FineH,Z_FineH,'k--','linewidth',3)                              
%                             
%         end
%     end                    
    figureCount = figureCount+1;                    
end      
