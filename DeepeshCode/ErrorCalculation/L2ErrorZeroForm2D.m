function errorL2 = L2ErrorZeroForm2D(zeroFormDiscrete, funct, phi, g, gridType)

% L2 Error Norm for zero forms.
% errorL2 = L2ErrorZeroForm2D(zeroFormDiscrete, funct, phi, g, gridType)
% L2ErrorNorm = (<error, error>)^{1/2}
% where 
% error = (computed - exact)
%
% Inputs:
%
%     phi                   :: mappings from reference domain to physical
%     funct                 :: zero form function
%     g                     :: square-root of metric tensor determinant
%     gridType              :: Mesh type
%     zeroFormDiscrete      :: reduced zero forms [nElements X (p+1)^2]
%
% Outputs:
%
%     errorL2                :: error in reconstruction
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $2/03/2012$

    % number of elements
    nElements = size(g,1);
    
    % order of mesh
    p = sqrt(size(zeroFormDiscrete,1))-1;
    
    % quadrature order
    pint = p+3;

    % Quadrature Nodes and weights
    [quadNodes, quadWeights] = GaussQuad(pint);
    QuadWeights = kron(quadWeights,quadWeights);

    % Reconstruct zeroFormDiscrete at quadNodes
    [zeroFormReconstructed meshReconstructionXi meshReconstructionEta] = ReconstructZeroForms2D(quadNodes, quadNodes, phi, zeroFormDiscrete, gridType);

    errorL2 = 0;

    % Reduce f at quadNodes
    for element = 1:nElements

        % Reduce f at quadNodes
        [meshReconstructionX meshReconstructionY] = phi{element}(meshReconstructionXi,meshReconstructionEta);
        zeroFormExactElement = funct(meshReconstructionX,meshReconstructionY);

        % Metric Tensor
        evaluatedG = g{element}(meshReconstructionXi,meshReconstructionEta);

        % Integrate square of difference for each element
        fDiffSqIntegral = QuadWeights'*(((zeroFormExactElement(:)-zeroFormReconstructed(:,element)).^2).*evaluatedG(:));

        errorL2 = errorL2 + fDiffSqIntegral;

    end

    %L2 error norm
    errorL2 = sqrt(errorL2);
    
end