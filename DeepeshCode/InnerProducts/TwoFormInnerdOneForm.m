function innerProd = TwoFormInnerdOneForm(p,g,nodesType,sparseFlag,varargin)
%TwoFormInnerdOneForm Computes all the inner products of the 2-form
%                    basis functions with the exterior derivative of 1-form
%                    basis functions
%
%   innerProd = TwoFormInnerdOneForm(p,g,nodesType,sparseFlag,~quadOrder,~quadType)
%
%   Where:
%       p        :: the order of the basis functions
%       g        :: the square root of the determinant of the g_{ij} metric,
%                   that is, the Jacobian of the mapping
%       nodesType  :: defines the type of nodes to use (Lobatto or EGauss)
%       sparseFlag :: defines if output should be sparse or not
%
%   Optional
%       quadOrder :: defines the order of integration to use (the same for
%                    both xi and eta directions)
%       quadType  :: defines the type of integration to use (the same for
%                    both xi and eta directions)
%
%   Returns the inner products:
%       < \omega_{i}, dxi_{n} > and < \omega_{i}, deta_{n} >
%
%   Where the \omega_{i} are 2-form basis functions and dxi_{n} and
%   deta_{n} are 1-form basis functions.

%   Copyright 2011 Artur Palha
%   $ Revision: 2.0 $  $ Date: 2011/11/01 $

    % check if nodesType is a valid one
    if ~TestPolyType(nodesType)
        disp(sprintf(':: %s :: is not a valid type of nodes', nodesType));
        return
    end

    % check for the number of inputs and define the inputs if default ones
    % are to be used
    if nargin > 4
        quadOrder = varargin{1};
        quadType = varargin{2};
        
        % check if quadType is a valid one
        if ~TestPolyType(quadType)
            disp(sprintf(':: %s :: is not a valid type of quadrature', quadType));
            return
        end
    else
        % default inputs
        quadOrder = p;
        quadType = 'Gauss';
    end
    
    %% Compute the inner products of 2-forms
    twoFormInnerProduct = TwoFormInnerTwoForm(p,g,nodesType,sparseFlag,quadOrder,quadType);
    
    %% Compute the exterior derivative of 1-forms
    dOneForm = dOne(p);
    
    %% Compute the inner product
    innerProd = twoFormInnerProduct * dOneForm;
end