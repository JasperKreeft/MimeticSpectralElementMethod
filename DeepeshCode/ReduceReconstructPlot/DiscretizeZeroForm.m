function zeroFormDiscrete = DiscretizeZeroForm(f, phi, p, gridType, varargin)

% Discretize zero-forms.
%
%   zeroFormDiscrete = DiscretizeZeroForm(f, phi, p, gridType)
%
%   Where:
%
%       p                 :: order of mesh
%       phi               :: mappings from reference domain to physical
%       f                 :: analytical zero-form function
%       gridType          :: the type of grid
%
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.0 $  $Date: 4/3/2012 $

    if size(varargin,2)
        time = varargin{1};
    else
        time = [];
    end

    % number of elements
    nElements = size(phi,1);

    % Zero Form Points in 1D
    xi = eval(sprintf('%sQuad(%s)', strtrim(gridType), 'p'));
    eta = xi;
    
    % Zero Form Points in 2D
    [Xi Eta] = meshgrid(xi,eta);
    
    % allocate memory space for the result of the discretization
    zeroFormDiscrete = zeros((p+1)*(p+1),nElements);
    
    if (size(time,1))
        for element = 1:nElements

            [x y] = phi{element}(Xi,Eta);

            tempZeroForm = f(x,y,time);

            zeroFormDiscrete(:,element) = tempZeroForm(:);

        end
    else
        for element = 1:nElements

            [x y] = phi{element}(Xi,Eta);

            tempZeroForm = f(x,y);

            zeroFormDiscrete(:,element) = tempZeroForm(:);

        end
    end
    
    clear tempZeroForm
    
end