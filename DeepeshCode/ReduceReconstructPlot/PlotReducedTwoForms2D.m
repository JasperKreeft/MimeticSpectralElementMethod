function PlotReducedTwoForms2D(twoFormDiscrete,phi,p,gridType,figureTwoFormDiscrete)

    nElements = size(phi,1);
    
    % GridNodes
    gridNodes = eval(sprintf('%sQuad(%s)', strtrim(gridType), 'p'));

    % Prepare twoFormBoundaries
    gridXi = rectpulse([gridNodes(1:p) gridNodes(2:p+1)],p);
    gridEta = repmat([gridNodes(1:p) gridNodes(2:p+1)],p,1);
    
    
    for element = 1:nElements
       
        [meshX meshY] = phi{element}(gridXi,gridEta);
        
        figure(figureTwoFormDiscrete)
        for twoForm = 1:p*p
            
            [x y] = meshgrid(meshX(twoForm,:),meshY(twoForm,:));
            z = twoFormDiscrete(twoForm,element)*ones(size(x));
            
            surf(x,y,z)
            hold on
            
        end
        grid on
        
    end


end