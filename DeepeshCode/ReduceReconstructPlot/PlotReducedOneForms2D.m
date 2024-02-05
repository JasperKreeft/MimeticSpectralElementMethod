function PlotReducedOneForms2D(oneFormDiscrete,phi,p,gridType,figureOneFormReduced)

    nElements = size(phi,1);
    zMax = max(max(oneFormDiscrete));
    zMin = min(min(oneFormDiscrete));
    
    % Map [zMin,zMax] -> [-1,1]
    % z = c*phi + d
    c = 0.5*(zMax - zMin);
    d = 0.5*(zMax + zMin);
    
    
    % Normalized zData
    normalizedZData = (oneFormDiscrete - d)/c;
    
    % GridNodes
    gridNodes = eval(sprintf('%sQuad(%s)', strtrim(gridType), 'p'));

    % Prepare OneFormBoundaries
    % Xi componenets
    gridXiXi = rectpulse([gridNodes(1:p,1) gridNodes(2:p+1,1)],p+1);
    gridEtaXi = repmat(repmat(gridNodes,1,2),p,1);
    if (size(gridXiXi)~=size(gridEtaXi)) % p = 1
        gridXiXi = reshape(gridXiXi,size(gridEtaXi));
    end
    % Eta components
    gridXiEta = rectpulse(repmat(gridNodes,1,2),p);
    gridEtaEta = repmat([gridNodes(1:p) gridNodes(2:p+1)],p+1,1);
    
    for element = 1:nElements
        
        [meshXXi meshYXi] = phi{element}(gridXiXi,gridEtaXi);
        [meshXEta meshYEta] = phi{element}(gridXiEta,gridEtaEta);
        
        figure(figureOneFormReduced)
        for oneFormXi = 1:p*(p+1)
            
            x = meshXXi(oneFormXi,:)';
            y = meshYXi(oneFormXi,:)';
            z = oneFormDiscrete(oneFormXi,element)*ones(2,1);
            
            weightRGB = LobattoPoly(normalizedZData(oneFormXi,element),2);
            color = fliplr((max(weightRGB,0))');
            
            plot3(x, y, z,'-','Color',color);
            hold on
            
        end
        
        for oneFormEta = p*(p+1)+1:2*p*(p+1)
            
            x = meshXEta(oneFormEta-p*(p+1),:)';
            y = meshYEta(oneFormEta-p*(p+1),:)';
            z = oneFormDiscrete(oneFormEta,element)*ones(2,1);
            
            weightRGB = LobattoPoly(normalizedZData(oneFormEta,element),2);
            color = fliplr((max(weightRGB,0))');
            
            plot3(x, y, z,'-','Color',color);
            hold on
            
        end
        grid on
        
        
    end

end