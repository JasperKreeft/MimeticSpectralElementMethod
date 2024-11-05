%% Plotten %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('PHI','var')
    plotPotential
end

if exist('Q','var')
    plotFlux(Meshp,qx,qy,qMag,qx_interp,qy_interp,q_interp, ...
             FunctionType,Domain)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%