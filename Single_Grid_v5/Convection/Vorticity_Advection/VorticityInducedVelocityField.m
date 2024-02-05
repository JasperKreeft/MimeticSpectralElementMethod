function [U,V,PSI] = VorticityInducedVelocityField(Wrhs,M0,PM_L,PM_U,NG)

disp('Solve vorticity induced velocity field')

global globalnr_0 globalnr_1v globalnr_1h

% Curl of stream-function

RHS = M0*Wrhs;

y   = PM_L\RHS;
PSI = PM_U\y;

% Gradient values
UV = NG*PSI;

U   = UV(globalnr_1v);
V   = UV(globalnr_1h);

if nargout==3
    PSI = PSI(globalnr_0);
end