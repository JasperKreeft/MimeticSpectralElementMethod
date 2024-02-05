function basis = GaussPoly2D(xRecPara,p,gridType)

% xRecPara = linspace(-1,1,10);
% p = 3;

nP = length(xRecPara);

ZeroFormBasis = GaussPoly(xRecPara,p);

basis = kron(phiX)