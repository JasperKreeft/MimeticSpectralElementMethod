clear all
close all
clc

%%% ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo %%%
%%% multi-element test case for 0-forms on the standard domain in 2D    %%%
%%% includes bc and mapping                                             %%%
%%% ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo %%%

%-------------------------------------------------------------------------%
% settings                                                                %
%-------------------------------------------------------------------------%

hxOrder = 2;      % element order x-dir.
hyOrder = 1;      % element order y-dir.
setEqualH = 1;    % equal order in x- en y-dir. > 1 = equal, 0 = not equal
% hxOrder = 1:1:3;    % element order x-dir. (max=7) 
% hyOrder = 1:1:2;    % element order y-dir. (max=7)
% setEqualH = 0;      % equal order in x- en y-dir. > 1 = equal, 0 = not equal        

pxOrder = 6;      % element order x-dir.
pyOrder = 6;      % element order y-dir.
setEqualP = 1;    % equal order in x- en y-dir. > 1 = equal, 0 = not equal
% pxOrder = [2 4 6 8];% element order x-dir. (max=7)  
% pyOrder = [2 4 6 8];% element order y-dir. (max=7) 
% pxOrder = 1:7;% element order x-dir.  
% pyOrder = [2 4 6 8];% element order y-dir.
% setEqualP = 1;      % equal order in x- en y-dir. > 1 = equal, 0 = not equal

plotNumberingFlag = 1;	% numbering > 1 = plot, 0 = no plot
plotResultsFlag = 1;   	% results > 1 = plot, 0 = no plot
plotErrorFlag = 0;     	% error > 1 = plot, 0 = no plot
plotCondFlagP = 0;      % conditioning > 1 = plot, 0 = no plot
plotCondFlagH = 0;      % conditioning > 1 = plot, 0 = no plot 

printFlag = 1;	% 1 = print, 0 = no print
saveFlag = 0;	% 1 = save, 0 = no save

mapMethod = 'scale';	% mapping method (NOTE: only use 'scale')
mapVar = [-1 0 -1 1];  	% mapping variables; [xLeft xRight yLeft yRight]

funcPar = [1 1 1];	% [function, function parameters]
fineGridN = 100;   	% refinement level of fine grid
errorQuad = 100;   	% integration order norm quadrature

BCs = [1 1 1 1];	% [S(1) N(2) W(3) E(4)], 1 = dirichlet, 2 = neumann

%-------------------------------------------------------------------------%
% variable initialization                                                 %
%-------------------------------------------------------------------------%

figureCount = 1;        % figure counter

nPx = length(pxOrder);      % length of order vector in x-dir.
if setEqualP == 1
    nPy = 1;                % length of order vector in y-dir.
elseif setEqualP == 0
    nPy = length(pyOrder);  % length of order vector in y-dir.
end
nPxy = nPx*nPy;             % total number of order combinations

nHx = length(hxOrder);      % length of element vector in x-dir.
if setEqualH == 1
    nHy = 1;                % length of element vector in y-dir.
elseif setEqualH == 0
    nHy = length(hyOrder);  % length of element vector in y-dir.
end
nHxy = nHx*nHy;             % total number of element combinations

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

xBounds = mapVar(1:2);	% since the domain is only scaled
yBounds = mapVar(3:4);	% since the domain is only scaled
xElementCoords = xBounds;   
yElementCoords = yBounds;  
if plotNumberingFlag == 1
    PlotGlobalNumberingZeroFormTwoD([hxOrder hyOrder],[pxOrder pyOrder],xBounds,yBounds,1)
    PlotGlobalNumberingOneFormTwoD([hxOrder hyOrder],[pxOrder pyOrder],xBounds,yBounds,[1 1])
    PlotGlobalNumberingTwoFormTwoD([hxOrder hyOrder],[pxOrder pyOrder],xBounds,yBounds,1)
    figureCount = figureCount+3;
end

%-------------------------------------------------------------------------%
% calculate exact solution (fine grid)                                    %
%-------------------------------------------------------------------------%

fineGridRef = linspace(-1,1,fineGridN);                                     % finegrid    
[xFine yFine] = mappingMethod(mapMethod,mapVar,fineGridRef,fineGridRef);	% interpolation mesh x-dir./y-dir.
U_Fine = funcTwoD(xFine,yFine,0,funcPar,1);                                 % exact solution (matrix)

%-------------------------------------------------------------------------%
% allocate storage space                                                  %
%-------------------------------------------------------------------------%

if setEqualH == 1
    error = zeros(length(hxOrder),1,length(pxOrder),length(pyOrder),2);
elseif setEqualP == 1
    error = zeros(length(hxOrder),length(hyOrder),length(pxOrder),1,2);
elseif setEqualH == 1 && setEqualP == 1
    error = zeros(length(hxOrder),1,length(pxOrder),1,2);    
elseif setEqualH == 0 && setEqualP == 0
    error = zeros(length(hxOrder),length(hyOrder),length(pxOrder),length(pyOrder),2);
end   

%-------------------------------------------------------------------------%
% set up and solve systems                                                %
%-------------------------------------------------------------------------%

tAll = tic;
elementCountX=1;
for hx=hxOrder
    elementCountY=1;
    
    %%% set equal order
    if setEqualH == 1
        hyOrder = hx;
    end     
    
    for hy=hyOrder
    
        xElementCoordsDefault = -1:2/hx:1;
        yElementCoordsDefault = -1:2/hy:1;
        [xElementCoords yElementCoords] = mappingMethod(mapMethod,mapVar,xElementCoordsDefault,yElementCoordsDefault);
    
        orderCountX=1;
        for px=pxOrder  
            orderCountY=1; 
            
         	%%% set equal order
            if setEqualP == 1
                pyOrder = px;
            end 
            
            for py=pyOrder  
                tLoop = tic;

                %---------------------------------------------------------%
                % standard element mesh and polynomials                   %
                %---------------------------------------------------------%

                [refMeshNodesX refMeshWeightsX] = GLLnodes(px);	% mesh reference nodes and weights (x-dir.)
                [refMeshNodesY refMeshWeightsY] = GLLnodes(py);	% mesh reference nodes and weights (y-dir.) 
                [hnx enx] = MimeticpolyVal(refMeshNodesX,px,1);	% coarse mesh interpolants (x-dir.) 
                [hny eny] = MimeticpolyVal(refMeshNodesY,py,1);	% coarse mesh interpolants (y-dir.)   
                [hfx efx] = MimeticpolyVal(fineGridRef,px,1);	% fine mesh interpolants (x-dir.)
                [hfy efy] = MimeticpolyVal(fineGridRef,py,1);	% fine mesh interpolants (y-dir.)

                %---------------------------------------------------------%
                % mesh geometry                                           %
                %---------------------------------------------------------%

                [X_MeshZeroForms Y_MeshZeroForms] = ZeroFormCoordsTwoD([hx hy],[px py],xBounds,yBounds); % 0-form coordinates
                [X_MeshOneForms Y_MeshOneForms] = OneFormCoordsTwoD([hx hy],[px py],xBounds,yBounds);    % 1-form coordinates
                [X_MeshTwoForms Y_MeshTwoForms] = TwoFormCoordsTwoD([hx hy],[px py],xBounds,yBounds);    % 2-form coordinates             
                
                %---------------------------------------------------------%
                % mesh topology                                           %
                %---------------------------------------------------------%    

                D10x = dZeroXiTwoD([px py]);	% incidence matrix x-dir.
                D10y = dZeroEtaTwoD([px py]);	% incidence matrix y-dir.

                ZeroForms = GlobalNumberingZeroFormTwoD([hx hy],[px py]);   % 0-forms numbering
                OneForms  = GlobalNumberingOneFormTwoD([hx hy],[px py]);    % 1-forms numbering
                TwoForms  = GlobalNumberingTwoFormTwoD([hx hy],[px py]);  	% 2-forms numbering    

                nNodes    = max(max(ZeroForms));	% # of 0-forms
                nLines    = max(max(OneForms));     % # of 1-forms
                nSurfaces = max(max(TwoForms));     % # of 2-forms

                %---------------------------------------------------------%
                % storage                                                 %
                %---------------------------------------------------------%

                A = zeros(nNodes,nNodes);	% system matrix
                b = zeros(nNodes,1);        % system vector

                boundaryNodesGn = cell(4,1);	% global numbering boundary nodes
                boundaryValues = cell(4,1);     % boundary values
                boundaryNodesLn = cell(4,1);   	% local numbering boundary nodes   
                boundaryElements = cell(4,1);   % global numbering boundary elements    

                plotExactStr = cell(2,1);	% exact plot title
                plotNodalStr = cell(4,1);	% nodal plot title
                plotInterpStr = cell(4,1);	% interpolation plot title

                %---------------------------------------------------------%
                % matrix system build-up                                  %
                %---------------------------------------------------------%

                for elementX=1:hx
                    for elementY=1:hy
                        
                        %%% global element numbering
                        element = elementX+hx*(elementY-1);
                        
                        %%% physical coordinates
                        xMesh = X_MeshZeroForms(element,:);
                        yMesh = Y_MeshZeroForms(element,:);

                        %%% rhs
                        fMesh = -funcTwoD(xMesh,yMesh,2,funcPar,3);

                        %%% jacobian
                        mapVarEl = [xElementCoords(elementX) xElementCoords(elementX+1) ...
                                    yElementCoords(elementY) yElementCoords(elementY+1)];                                
                        [J D] = jacobianTwoD(mapMethod,mapVarEl);   % NOTE: J is just a constant for 'scale'
                        dx_dxi = D(1,1);                           
                        dy_deta = D(2,2);   

                        %%% building blocks
                        Wx = spdiags(refMeshWeightsX,0,px+1,px+1);  % integration weights x-dir.
                        Wy = spdiags(refMeshWeightsY,0,py+1,py+1);  % integration weights y-dir.
                        M0x = hnx*Wx*hnx';                         	% lagrange integration x-dir.
                        M0y = hny*Wy*hny';                        	% lagrange integration y-dir.
                        M1x = enx*Wx*enx';                        	% edge integration x-dir.
                        M1y = eny*Wy*eny';                        	% edge integration y-dir.

                        A1 = D10x'*kron(M0y,M1x)*D10x/J*dy_deta^2;  % contribution in x-dir.
                        A2 = D10y'*kron(M0x,M1y)*D10y/J*dx_dxi^2;   % contribution in y-dir.
                        RHS = kron(M0y,M0x)*fMesh*J;                % right hand side

                        A(ZeroForms(element,:),ZeroForms(element,:)) = ...
                            A(ZeroForms(element,:),ZeroForms(element,:))+A1+A2;   % fill system matrix               
                        b(ZeroForms(element,:),1) = b(ZeroForms(element,:),1)+RHS;                        % fill rhs vector 
                    
                    end
                end

                %---------------------------------------------------------%
                % numbering of boundary nodes/elements                    %
                %---------------------------------------------------------%
                
                %%% local numbering 
                boundaryNodesLn{1} = 1:(px+1);                          % lower boundary
                boundaryNodesLn{2} = boundaryNodesLn{1}+(py*(px+1));	% upper boundary
                boundaryNodesLn{3} = 1:(px+1):(1+py*(px+1));            % left boundary
                boundaryNodesLn{4} = boundaryNodesLn{3}+px;             % right boundary

                %%% boundary elements
                boundaryElements{1} = 1:hx;                             % lower boundary          
                boundaryElements{2} = boundaryElements{1}+(hy-1)*hx;    % upper boundary
                boundaryElements{3} = 1:hx:(1+hx*(hy-1));               % left boundary
                boundaryElements{4} = boundaryElements{3}+(hx-1);      	% right boundary          
                
                %%% global numbering 
                boundaryNodesGn{1} = ZeroForms(boundaryElements{1},boundaryNodesLn{1})'; 	% lower boundary 
                boundaryNodesGn{2} = ZeroForms(boundaryElements{2},boundaryNodesLn{2})';	% upper boundary 
                boundaryNodesGn{3} = ZeroForms(boundaryElements{3},boundaryNodesLn{3})';	% left boundary
                boundaryNodesGn{4} = ZeroForms(boundaryElements{4},boundaryNodesLn{4})';	% right boundary              
                
                %---------------------------------------------------------%
                % apply boundary conditions 	                          %
                %---------------------------------------------------------%        

                %%% dirichlet
                for nD = dirichletSides           
                    boundaryValues{nD} = ...
                        funcTwoD(X_MeshZeroForms(boundaryElements{nD},boundaryNodesLn{nD}),...
                        Y_MeshZeroForms(boundaryElements{nD},boundaryNodesLn{nD}),0,funcPar,3);
                end
                %%% neumann 
                for nN = neumannSides
                    ...
                end
                %%% call bc function
                [A b nodesOut firstVal dirichletValues] = ...
                    bcZeroFormTwoD(BCs,A,b,boundaryNodesGn,boundaryValues);

                %---------------------------------------------------------%
                % solve                                                   %
                %---------------------------------------------------------%

                error(elementCountX,elementCountY,orderCountX,orderCountY,2) = cond(A);
                A = sparse(A);
                b = sparse(b);
                x = A\b;

                %---------------------------------------------------------%
                % reconstruct solution vector                             %
                %---------------------------------------------------------%    

                if dirichletFlag
                    nodesIn = 1:nNodes; nodesIn(nodesOut) = [];
                    sol = zeros(nNodes,1); 
                    sol(nodesIn) = x;
                    sol(nodesOut) = dirichletValues(firstVal);        
                else
                    sol = x;
                end

                %---------------------------------------------------------%
                % error determination                                     %
                %---------------------------------------------------------%    
                
                for elementX=1:hx
                    for elementY=1:hy  
                        
                        element = elementX+hx*(elementY-1);           
   
                        mapVarEl = [xElementCoords(elementX) xElementCoords(elementX+1) ...
                                    yElementCoords(elementY) yElementCoords(elementY+1)];

                        elementSol = sol(ZeroForms(element,:));
                        elementSol = reshape(elementSol,px+1,py+1)';   
                        
                        error(elementCountX,elementCountY,orderCountX,orderCountY,1) = ...
                            error(elementCountX,elementCountY,orderCountX,orderCountY,1)+...
                            normsTwoD(elementSol,funcPar,0,mapMethod,mapVarEl,errorQuad);     
                    end
                end

                %---------------------------------------------------------%
                % print output                                            %
                %---------------------------------------------------------%        

                time = toc(tLoop);
                if printFlag == 1
                    fm = '%s %2.0f \t %s %2.0f \t %s %2.0f \t %s %2.0f \t %s %1.3e \t %s %1.3e \t %s %1.4f \n';
                    fprintf(fm,'hx =',hx,'hy =',hy,'px =',px,'py =',py,'L2 =',...
                        error(elementCountX,elementCountY,orderCountX,orderCountY,1),...
                        'cond =',error(elementCountX,elementCountY,orderCountX,orderCountY,1),'t =',time);
                end                 
                
                %---------------------------------------------------------%
                % plot results                                            %
                %---------------------------------------------------------%        
%%
                if plotResultsFlag == 1
                    
                    %%% plotting the exact solution    
                    figure(figureCount),subplot(1,3,1), hold on                                                           
                    surf(xFine,yFine,U_Fine)
                    shading flat, axis tight, grid on, axis equal,
                    plotExactStr{1} = 'exact';
                    plotExactStr{2} = strcat('fineGridN=',num2str(fineGridN));
                    title(plotExactStr), xlabel('x'), ylabel('y'), zlabel('z'), view(3)

                    %%% loop over elements, interpolate and plot results
                    for elementX=1:hx
                        for elementY=1:hy
                            
                            %%% global element numbering
                            element = elementX+hx*(elementY-1);
                            
                            %%% plot element corner nodes
                            subplot(1,3,2), hold on        
                            xEl = [xElementCoords(elementX) xElementCoords(elementX+1)];
                            yEl = [yElementCoords(elementY) yElementCoords(elementY+1)];
                            Z_ElementValues = funcTwoD(xEl,yEl,0,funcPar,1);
                            [X_ElementCoords Y_ElementCoords] = meshgrid(xEl,yEl);
                            plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2), view(3)                                                         
                            
                            %%% construct solution matrix
                            elementSol = sol(ZeroForms(element,:));
                            elementSol = reshape(elementSol,px+1,py+1)';
                            interpSol = hfy'*elementSol*hfx;            

                            %%% coloring
                            myColor = zeros(length(sol),1);
                            vmin = min(sol); vmax = max(sol); vrange = vmax-vmin';
                            for k=1:length(sol)
                                myColor(k) = (sol(k)-vmin)/vrange+1;
                            end

                            %%% nodal coordinates   
                            xScat = X_MeshZeroForms(element,:)';
                            yScat = Y_MeshZeroForms(element,:)';                            

                            %%% plot solution nodes
                            myscat = scatter3(xScat,yScat,sol(ZeroForms(element,:)),25,myColor(ZeroForms(element,:)),'filled');
                            plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2)
                            axis tight, grid on, axis equal, title('nodal'), view(3)
                            plotNodalStr{1} = 'nodal';
                            plotNodalStr{2} = strcat('h_x=',num2str(hx),', h_y=',num2str(hy));                            
                            plotNodalStr{3} = strcat('p_x=',num2str(px),', p_y=',num2str(py));
                            plotNodalStr{4} = strcat('BC=[',num2str(BCs),']');
                            title(plotNodalStr), xlabel('x'), ylabel('y'), zlabel('z')

                            %%% fine coordinates
                            mapVarEl = [xEl yEl];
                            [xFineEl yFineEl] = mappingMethod(mapMethod,mapVarEl,fineGridRef,fineGridRef); 
                            
                            %%% plot interpolated solution
                            subplot(1,3,3), hold on
                            surf(xFineEl,yFineEl,interpSol); 
                            plot3(X_ElementCoords,Y_ElementCoords,Z_ElementValues,'k+','linewidth',2)
                            shading flat, axis tight, grid on, axis equal, view(3)
                            plotInterpStr{1} = 'interpolated';
                            plotInterpStr{2} = strcat('h_x=',num2str(hx),', h_y=',num2str(hy));                            
                            plotInterpStr{3} = strcat('p_x=',num2str(px),', p_y=',num2str(py));
                            plotInterpStr{4} = strcat('BC=[',num2str(BCs),']');            
                            title(plotInterpStr), xlabel('x'), ylabel('y'), zlabel('z')
                       
                            %%% plot element boundaries
                            X_FineV = repmat(xElementCoords',fineGridN,1);
                            Y_FineV = repmat(yFine,1,hx+1);
                            Z_FineV = funcTwoD(X_FineV,Y_FineV,0,funcPar,3);
                            plot3(X_FineV,Y_FineV,Z_FineV,'k--','linewidth',3)

                            Y_FineH = repmat(yElementCoords',fineGridN,1);
                            X_FineH = repmat(xFine,1,hy+1);
                            Z_FineH = funcTwoD(X_FineH,Y_FineH,0,funcPar,3);
                            plot3(X_FineH,Y_FineH,Z_FineH,'k--','linewidth',3)                              

                            %%% cutting planes
%                             for i=2:hx
%                                 x = [xElementCoords(i);xElementCoords(i)];
%                                 y = [yElementCoords(1);yElementCoords(hy+1)];
%                                 [X Y] = meshgrid(x,y);
%                                 zMin = min(min(U_Fine));
%                                 zMax = max(max(U_Fine));
%                                 Z = repmat([zMin zMax],2,1);
%                                 surf(X,Y,Z)
%                                 alpha(.25)
%                             end
%                             for i=2:hy
%                                 x = [xElementCoords(1);xElementCoords(hx+1)];                                
%                                 y = [yElementCoords(i);yElementCoords(i)];
%                                 [X Y] = meshgrid(x,y);
%                                 zMin = min(min(U_Fine));
%                                 zMax = max(max(U_Fine));
%                                 Z = repmat([zMin; zMax],1,2);
%                                 surf(X,Y,Z)
%                                 alpha(1)
%                             end
                            
                        end
                    end                    
                    figureCount = figureCount+1;                    
                end      

                orderCountY = orderCountY+1;
            end       
            orderCountX = orderCountX+1;
        end
        elementCountY = elementCountY+1;
    end
    elementCountX = elementCountX+1;
end

tFinal = toc(tAll);
fm = '%s %1.4f \n';
fprintf(fm,'total run time =',tFinal);

%%
%-------------------------------------------------------------------------%
% plot error                                                              %
%-------------------------------------------------------------------------%

if plotErrorFlag == 1    
    %%% coloring data
    colorOpts = [0 0 1; % blue
                 0 1 0; % green
                 1 0 0; % red                         
                 0 1 1; % cyan
                 1 0 1; % magenta
                 1 1 0; % yellow 
                 0 0 0];% black          
    ind = [3 0;
           2 0;
           1 0;       
           2 0;
           1 0;
           1 0;
           1 1];    
    % 2nd row: 0=- 1=+
    % myMarkers = {'--^','--o','--*','--x','--s','--d','--+'};
    myMarkers = {'-^','-o','-s','-d','-*','-x','-+'};  
    
    %%% plotting options 
    figure(figureCount), hold on
    if setEqualH == 1 && setEqualP == 1 % equal h and p        
        errorPlotHndl = zeros(nHx,1);
        legendtext = {zeros(nHx,1)};         
        for elX=1:nHx
            legendtext(elX) = {strcat('h = [',num2str(hxOrder(elX)),',',...
            	num2str(hxOrder(elX)),']')};
            val = 1-.5/nHx*(elX-1);
            myBlue = [0 0 val];
            if elX==1
                errorPlotHndl(elX) = plot(pxOrder,squeeze(error(elX,1,:,1,1)),'o-','color',myBlue,'linewidth',2,'MarkerEdgeColor','g');
            elseif elX==nHx
                errorPlotHndl(elX) = plot(pxOrder,squeeze(error(elX,1,:,1,1)),'o-','color',myBlue,'linewidth',2,'MarkerEdgeColor','r');
            else
                errorPlotHndl(elX) = plot(pxOrder,squeeze(error(elX,1,:,1,1)),'o-','color',myBlue,'linewidth',2);
            end
        end       
        gridLegend(errorPlotHndl,nHx,legendtext,'location','southoutside','orientation','horizontal');
        % legend(legendtext,'location','EastOutside')
        xlabel('p')        
    elseif setEqualH == 1 % equal h
        errorPlotHndl = zeros(nHx*nPx,1);
        legendtext = {zeros(nHx*nPx,1)};
        for elX=1:nHx
            for i = 1:nPx                
                plotCount = (elX-1)*nPx+i;
                legendtext(plotCount) = {strcat('h = [',num2str(hxOrder(elX)),',',...
                	num2str(hxOrder(elX)),']',', p_x = ',num2str(pxOrder(i)))};
                myColor = colorOpts(elX,:);
                if ind(elX,2) == 0
                    val = 1-.5/nPx*(i-1);
                elseif ind(elX,2) == 1
                    val = .5/nPx*(i-1);
                end                
                myColor(ind(elX,1)) = val;    
                errorPlotHndl(plotCount)=plot(pyOrder,squeeze(error(elX,1,i,:,1)),myMarkers{elX},'color',myColor,'linewidth',2);
            end
        end
        gridLegend(errorPlotHndl,nPx,legendtext,'location','southoutside','orientation','horizontal');
        % legend(legendtext,'location','EastOutside') 
        xlabel('p_y')
    elseif setEqualP == 1 % equal p
        errorPlotHndl = zeros(nHxy,1);
        legendtext = {zeros(nHxy,1)};        
        for elX=1:nHx
            for elY=1:nHy        
                plotCount = elY+nHy*(elX-1);
                legendtext(plotCount) = {strcat('h = [',num2str(hxOrder(elX)),',',...
                	num2str(hyOrder(elY)),']')};
                myColor = colorOpts(plotCount,:); 
                errorPlotHndl(plotCount)=plot(pxOrder,squeeze(error(elX,elY,:,1,1)),myMarkers{plotCount},'color',myColor,'linewidth',2);
            end
        end  
        gridLegend(errorPlotHndl,nHxy,legendtext,'location','southoutside','orientation','horizontal');
        % legend(legendtext,'location','EastOutside')  
        xlabel('p')
    elseif setEqualH == 0 && setEqualP == 0 % general 
        errorPlotHndl = zeros(nHxy*nPx,1);
        legendtext = {zeros(nHxy*nPx,1)}; 
        for elX=1:nHx
            for elY=1:nHy
                hNumber = elY+nHy*(elX-1);
                for i = 1:nPx                
                    plotCount = (hNumber-1)*nPx+i;
                    legendtext(plotCount) = {strcat('h = [',num2str(hxOrder(elX)),',',...
                    	num2str(hyOrder(elY)),']',', p_x = ',num2str(pxOrder(i)))};
                    myColor = colorOpts(hNumber,:);
                    if ind(hNumber,2) == 0
                        val = 1-.5/nPx*(i-1);
                    elseif ind(hNumber,2) == 1
                        val = .5/nPx*(i-1);
                    end                
                    myColor(ind(hNumber,1)) = val;    
                    errorPlotHndl(plotCount)=plot(pyOrder,squeeze(error(elX,elY,i,:,1)),myMarkers{hNumber},'color',myColor,'linewidth',2);
                end
            end
        end   
        gridLegend(errorPlotHndl,nPx,legendtext,'location','southoutside','orientation','horizontal');
        % legend(legendtext,'location','EastOutside')  
        xlabel('p_y')
    end
    grid on, axis square, axis tight, title('L^2'), set(gca,'YScale','log')    
    figureCount=figureCount+1;
end  

%%
%-------------------------------------------------------------------------%
% plot cond (function of p)                                               %
%-------------------------------------------------------------------------%

if plotErrorFlag == 1    
    %%% coloring data
    colorOpts = [0 0 1; % blue
                 0 1 0; % green
                 1 0 0; % red                         
                 0 1 1; % cyan
                 1 0 1; % magenta
                 1 1 0; % yellow 
                 0 0 0];% black          
    ind = [3 0;
           2 0;
           1 0;       
           2 0;
           1 0;
           1 0;
           1 1];    
    % 2nd row: 0=- 1=+
    % myMarkers = {'--^','--o','--*','--x','--s','--d','--+'};
    myMarkers = {'-^','-o','-s','-d','-*','-x','-+'};  
    
    %%% plotting options 
    figure(figureCount), hold on
    if setEqualH == 1 && setEqualP == 1 % equal h and p        
        errorPlotHndl = zeros(nHx,1);
        legendtext = {zeros(nHx,1)};         
        for elX=1:nHx
            legendtext(elX) = {strcat('h = [',num2str(hxOrder(elX)),',',...
            	num2str(hxOrder(elX)),']')};
            val = 1-.5/nHx*(elX-1);
            myBlue = [0 0 val];
            if elX==1
                errorPlotHndl(elX) = plot(pxOrder,squeeze(error(elX,1,:,1,2)),'o-','color',myBlue,'linewidth',2,'MarkerEdgeColor','g');
            elseif elX==nHx
                errorPlotHndl(elX) = plot(pxOrder,squeeze(error(elX,1,:,1,2)),'o-','color',myBlue,'linewidth',2,'MarkerEdgeColor','r');
            else
                errorPlotHndl(elX) = plot(pxOrder,squeeze(error(elX,1,:,1,2)),'o-','color',myBlue,'linewidth',2);
            end
        end       
        gridLegend(errorPlotHndl,nHx,legendtext,'location','southoutside','orientation','horizontal');
        % legend(legendtext,'location','EastOutside')
        xlabel('p')        
    elseif setEqualH == 1 % equal h
        errorPlotHndl = zeros(nHx*nPx,1);
        legendtext = {zeros(nHx*nPx,1)};
        for elX=1:nHx
            for i = 1:nPx                
                plotCount = (elX-1)*nPx+i;
                legendtext(plotCount) = {strcat('h = [',num2str(hxOrder(elX)),',',...
                	num2str(hxOrder(elX)),']',', p_x = ',num2str(pxOrder(i)))};
                myColor = colorOpts(elX,:);
                if ind(elX,2) == 0
                    val = 1-.5/nPx*(i-1);
                elseif ind(elX,2) == 1
                    val = .5/nPx*(i-1);
                end                
                myColor(ind(elX,1)) = val;    
                errorPlotHndl(plotCount)=plot(pyOrder,squeeze(error(elX,1,i,:,2)),myMarkers{elX},'color',myColor,'linewidth',2);
            end
        end
        gridLegend(errorPlotHndl,nPx,legendtext,'location','southoutside','orientation','horizontal');
        % legend(legendtext,'location','EastOutside') 
        xlabel('p_y')
    elseif setEqualP == 1 % equal p
        errorPlotHndl = zeros(nHxy,1);
        legendtext = {zeros(nHxy,1)};        
        for elX=1:nHx
            for elY=1:nHy        
                plotCount = elY+nHy*(elX-1);
                legendtext(plotCount) = {strcat('h = [',num2str(hxOrder(elX)),',',...
                	num2str(hyOrder(elY)),']')};
                myColor = colorOpts(plotCount,:); 
                errorPlotHndl(plotCount)=plot(pxOrder,squeeze(error(elX,elY,:,1,2)),myMarkers{plotCount},'color',myColor,'linewidth',2);
            end
        end  
        gridLegend(errorPlotHndl,nHxy,legendtext,'location','southoutside','orientation','horizontal');
        % legend(legendtext,'location','EastOutside')  
        xlabel('p')
    elseif setEqualH == 0 && setEqualP == 0 % general 
        errorPlotHndl = zeros(nHxy*nPx,1);
        legendtext = {zeros(nHxy*nPx,1)}; 
        for elX=1:nHx
            for elY=1:nHy
                hNumber = elY+nHy*(elX-1);
                for i = 1:nPx                
                    plotCount = (hNumber-1)*nPx+i;
                    legendtext(plotCount) = {strcat('h = [',num2str(hxOrder(elX)),',',...
                    	num2str(hyOrder(elY)),']',', p_x = ',num2str(pxOrder(i)))};
                    myColor = colorOpts(hNumber,:);
                    if ind(hNumber,2) == 0
                        val = 1-.5/nPx*(i-1);
                    elseif ind(hNumber,2) == 1
                        val = .5/nPx*(i-1);
                    end                
                    myColor(ind(hNumber,1)) = val;    
                    errorPlotHndl(plotCount)=plot(pyOrder,squeeze(error(elX,elY,i,:,2)),myMarkers{hNumber},'color',myColor,'linewidth',2);
                end
            end
        end   
        gridLegend(errorPlotHndl,nPx,legendtext,'location','southoutside','orientation','horizontal');
        % legend(legendtext,'location','EastOutside')  
        xlabel('p_y')
    end
    grid on, axis square, axis tight, title('L^2'), set(gca,'YScale','log')    
    figureCount=figureCount+1;
end 

%%
%-------------------------------------------------------------------------%
% plot cond (function of h)                                               %
%-------------------------------------------------------------------------%

if plotCondFlagH == 1  
    %%% coloring data
    colorOpts = [0 0 1; % blue
                 0 1 0; % green
                 1 0 0; % red                         
                 0 1 1; % cyan
                 1 0 1; % magenta
                 1 1 0; % yellow 
                 0 0 0];% black   
    ind = [3 0;
           2 0;
           1 0;       
           2 0;
           1 0;
           1 0;
           1 1];   
     
    %%% plotting   
    figure(figureCount), hold on
    if setEqualH == 0 && setEqualP == 1  
        
        legendtext = {zeros(nPxy,1)};         
        conditionNumber = zeros(nHxy*nPxy,1);
        totalCount = 1;
        for elX=1:nHx
            for elY=1:nHy
                for i = 1:nPx                
                    conditionNumber(totalCount) = error(elX,elY,i,1,2);
                    totalCount = totalCount + 1;
                end
            end
        end  
        
        for i=1:nPxy                   
            surfData = conditionNumber(((1:nPxy:nHxy*nPxy)+(i-1)));            
            [X_Plot Y_Plot] = meshgrid(hxOrder,hyOrder);
            Z_Plot = reshape(surfData,nHy,nHx);  
            
            legendtext{i} = strcat('p = ',num2str(pxOrder(i)));
            
            myColorMatrix = zeros(nHy,nHx,3);           
            vmin = min(min(Z_Plot)); vmax = max(max(Z_Plot)); vrange = vmax-vmin;               
            for m = 1:nHy
                for n = 1:nHx
                    colorVal = colorOpts(i,:);
                    if ind(i,2) == 0
                    	colorVal(ind(i,1)) = 1-(Z_Plot(m,n)-vmin)/vrange/2;
                    else
                        colorVal(ind(i,1)) = (Z_Plot(m,n)-vmin)/vrange/2;
                        colorVal(1) = 0.5+(Z_Plot(m,n)-vmin)/vrange/2;
                        colorVal(2) = 0.5+(Z_Plot(m,n)-vmin)/vrange/2;
                        colorVal(3) = 0.5+(Z_Plot(m,n)-vmin)/vrange/2;
                    end
                    myColorMatrix(m,n,:) = colorVal;                    
                end
            end            
            surf(X_Plot,Y_Plot,Z_Plot,myColorMatrix);
            
        end
        grid on, view(3),axis square, axis tight, title('cond'), set(gca,'ZScale','log') 
        xlabel('hx'), ylabel('hy'), zlabel('cond')     
        legend(legendtext)
    end
    figureCount=figureCount+1;
end 

%%
%-------------------------------------------------------------------------%
% save data                                                               %
%-------------------------------------------------------------------------%
if saveFlag == 1
    mySettings.hx = hxOrder; mySettings.hy = hyOrder;
    mySettings.px = pxOrder; mySettings.py = pyOrder;
    mySettings.map = mapMethod; mySettings.mapvar = mapVar;
    mySettings.func = funcPar; mySettings.bc = BCs;
    mySettings.quad = errorQuad; mySettings.error = error;
    save('settings','mySettings');
end
