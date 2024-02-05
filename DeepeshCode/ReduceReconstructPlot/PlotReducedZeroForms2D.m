function PlotReducedZeroForms2D(zeroFormDiscrete,phi,gridType,figureZeroFormDiscrete)

    p = sqrt(size(zeroFormDiscrete,1))-1;

    % Zero Form Points
    xi = eval(sprintf('%sQuad(%s)', strtrim(gridType), 'p'));
    eta = xi;
    
    [Xi Eta] = meshgrid(xi,eta);

    nElements = size(phi,1);
    
    zMax = max(max(zeroFormDiscrete));
    zMin = min(min(zeroFormDiscrete));
    
    % Map [zMin,zMax] -> [-1,1]
    % z = c*phi + d
    c = 0.5*(zMax - zMin);
    d = 0.5*(zMax + zMin);
    
    % Normalized zData
    normalizedZData = (zeroFormDiscrete - d)/c;
    iMax = find(normalizedZData>1);
    normalizedZData(iMax) = ones(size(normalizedZData(iMax)));
    
    figure(figureZeroFormDiscrete)
    for element = 1:nElements
        
        [meshX meshY] = phi{element}(Xi,Eta);
        
        meshX = meshX(:);
        meshY = meshY(:);
        
        for zeroForm = 1:length(meshX)

            weightRGB = LobattoPoly(normalizedZData(zeroForm,element),2);
            color = fliplr((max(weightRGB,0))');
            
            plot3(meshX(zeroForm),meshY(zeroForm),zeroFormDiscrete(zeroForm,element),'.','Color',color);
            hold on
            
        end
        grid on
        
        
    end





end