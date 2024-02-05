function [PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,bc)

% int a^1 wedge star b^2     or     int a^{n-1} wedge star b^n

global N numElements numRows numColumns
global globalnr_1v globalnr_1h
global nr_1

if ~exist('numElements','var')
    numElements = 1;
end

Qbc    = [];
% boundary_points = [];
boundary_flux   = [];

PHIbc = zeros(nr_1,1);

for r=1:numRows
    for c=1:numColumns

if c==1
    ind1 = 1:N+1:N*(N+1);
    ind2 = c+(r-1)*numColumns;
    boundary_pointsL = globalnr_1v(ind1,ind2);
    if bc(1) == 1
        PHIbcL = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,'left','two');
        PHIbc(boundary_pointsL) = PHIbcL;
%         boundary_points = [ boundary_points ; boundary_pointsL ];
    else
        QbcL = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,'left','one');
        Qbc = [ Qbc ; QbcL' ];
        boundary_flux = [ boundary_flux ; boundary_pointsL ];
    end
end

if c==numColumns
    ind1 = N+1:N+1:N*(N+1);
    ind2 = r*numColumns;
    boundary_pointsR = globalnr_1v(ind1,ind2);
    if bc(2) == 1
        PHIbcR = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,'right','two');
        PHIbc(boundary_pointsR) = PHIbcR;
%         boundary_points = [ boundary_points ; boundary_pointsR ];
    else
        QbcR = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,'right','one');
        Qbc = [ Qbc ; QbcR' ];
        boundary_flux = [ boundary_flux ; boundary_pointsR ];
    end
end

if r==1
    ind1 = 1:N+1:N*(N+1);
    ind2 = c;
    boundary_pointsB = globalnr_1h(ind1,ind2);    
    if bc(3) == 1
        PHIbcB = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,'below','two');
        PHIbc(boundary_pointsB) = PHIbcB;
%         boundary_points = [ boundary_points ; boundary_pointsB ];
    else
        QbcB = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,'below','one');
        Qbc = [ Qbc ; QbcB ];
        boundary_flux = [ boundary_flux ; boundary_pointsB ];
    end
end

if r==numRows
    ind1 = N+1:N+1:N*(N+1);
    ind2 = (numRows-1)*numColumns+c;
    boundary_pointsA = globalnr_1h(ind1,ind2);
    if bc(4) == 1
        PHIbcA = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,'above','two');
        PHIbc(boundary_pointsA) = PHIbcA;
%         boundary_points = [ boundary_points ; boundary_pointsA ];
    else
        QbcA = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,'above','one');
        Qbc = [ Qbc ; QbcA ];
        boundary_flux = [ boundary_flux ; boundary_pointsA ];
    end
end

    end
end

interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];
