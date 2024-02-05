function [U,V,PSI] = VorticityInducedVelocityField_old(Wrhs,PSIbc,M0,Matrixrhs,PM_L,PM_U,NG,boundary_points,interior_points)

disp('Solve vorticity induced velocity field')

global globalnr_0 globalnr_1v globalnr_1h
global nr_0

% Curl of stream-function

RHS = M0*Wrhs;

if ~isempty(boundary_points)
    RHS = RHS - Matrixrhs*PSIbc;
end
RHS(boundary_points) = [];

y     = PM_L\RHS;
PSIin = PM_U\y;


PSI = zeros(nr_0,1);
PSI(interior_points) = PSIin;
PSI(boundary_points) = PSIbc;

% Gradient values
UV = NG*PSI;

U   = UV(globalnr_1v);
V   = UV(globalnr_1h);

if nargout==3
    PSI = PSI(globalnr_0);
end