phimin = -1; % min(phi);
phimax =  1; % max(phi);
umin    = -pi;
umax    =  pi;
velomin =   0;
velomax =  pi;

if exist('PHI','var')
plotPotential
end

if exist('U','var')
plotVelocity
end

if exist('Q','var')
plotFlux
end
