function [zeroFormReconstructed meshReconstructionXi meshReconstructionEta] = ReconstructZeroForms2D(discreteZeroForm, xiRecPara, etaRecPara, phi, gridType, varargin)

% Reconstruct zero forms.
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $2/03/2012$

    nElements = size(phi,1);
    discreteZeroForm = reshape(discreteZeroForm,size(discreteZeroForm,1)*size(discreteZeroForm,2)/nElements,nElements);    

    % gridType = 'Lobatto';
    p = sqrt(size(discreteZeroForm,1))-1;
    
    % Zero form bases
    ZeroFormBasis1DXi = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'xiRecPara', 'p'));
    ZeroFormBasis1DEta = eval(sprintf('%sPoly(%s,%s)', strtrim(gridType), 'etaRecPara', 'p'));

    % tensor product of 1D functions
    ZeroFormBasis2D = kron(ZeroFormBasis1DXi,ZeroFormBasis1DEta);

    % Reconstruct the function        
    zeroFormReconstructed = ZeroFormBasis2D'*discreteZeroForm;

    %Form Mesh
    [xi eta] = meshgrid(xiRecPara,etaRecPara);

    % mesh
    meshReconstructionXi = xi;
    meshReconstructionEta = eta;
    
end