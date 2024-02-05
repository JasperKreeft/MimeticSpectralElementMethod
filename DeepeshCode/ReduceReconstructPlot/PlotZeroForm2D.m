function PlotZeroForm2D(discreteZeroForm,phi,xBounds,yBounds,nReconstruction,gridType,figureNumber,varargin)
% Plot interpolations of zero-forms
%
%   PlotZeroForm2D(zeroFormDiscrete,phi,nReconstruction,gridType,figureZeroForm,xBounds,yBounds,varargin)
%
%   Where:
%
%       zeroFormDiscrete  :: the 0-form discretized
%       phi               :: a cell vector with the mappings from the
%                            parametric space to the physical space for each element
%       xBounds           :: the x bounds of the physical domain
%       yBounds           :: the y bounds of the physical domain
%       nReconstruction   :: the number of points to use in the x and y
%                            direction for the refinement
%       gridType          :: the type of grid
%       figureNumber      :: the number for the figure where to plot
%       varargin          :: 1. Contours (plot contours)
%                            2. n (number of contours)
%
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 4/3/2012 $


    xiRecPara = linspace(-1,1,nReconstruction);
    etaRecPara = xiRecPara;
    [zeroFormReconstructed meshReconstructionXi meshReconstructionEta] = ReconstructZeroForms2D(discreteZeroForm, xiRecPara, etaRecPara, phi, gridType);
    
    figure(figureNumber)
    axis([xBounds yBounds])
    hold on
    for element = 1:size(discreteZeroForm,2)
        
        [meshReconstructionX meshReconstructionY] = phi{element}(meshReconstructionXi,meshReconstructionEta);
        
        if (size(varargin,2))
            n = varargin{2};
            contour(meshReconstructionX,meshReconstructionY,reshape(zeroFormReconstructed(:,element),nReconstruction,nReconstruction),n)
            hold on
        else
            surf(meshReconstructionX,meshReconstructionY,reshape(zeroFormReconstructed(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            shading interp
            hold on
        end
    end
    grid on

end