function [PHIbc,Qbc,boundary_points_q,interior_points_q] = boundaryconditions(FunctionType,Domain,DomInfo,bc)

global globalnr_1v globalnr_1h
global nr_1

PHIbc_ = [];
boundary_points_u = [];
Qbc = [];
boundary_points_q = [];
if bc(1) == 1 
    PHIbcL = boundary_oneforms(FunctionType,Domain,DomInfo,'left','potential');
    PHIbc_ = [ PHIbc_ ; PHIbcL' ];
    boundary_points_u = [ boundary_points_u globalnr_1v(1,:) ];
else
    QbcL = boundary_oneforms(FunctionType,Domain,DomInfo,'left','flux');
    Qbc = [ Qbc ; QbcL' ];
    boundary_points_q = [ boundary_points_q globalnr_1v(1,:) ];
end
if bc(2) == 1
    PHIbcR = boundary_oneforms(FunctionType,Domain,DomInfo,'right','potential');
    PHIbc_ = [ PHIbc_ ; PHIbcR' ];
    boundary_points_u = [ boundary_points_u globalnr_1v(end,:) ];
else
    QbcR = boundary_oneforms(FunctionType,Domain,DomInfo,'right','flux');
    Qbc = [ Qbc ; QbcR' ];
    boundary_points_q = [ boundary_points_q globalnr_1v(end,:) ];
end
if bc(3) == 1
    PHIbcB = boundary_oneforms(FunctionType,Domain,DomInfo,'below','potential');
    PHIbc_ = [ PHIbc_ ; PHIbcB ];
    boundary_points_u = [ boundary_points_u globalnr_1h(:,1)' ];
else
    QbcB = boundary_oneforms(FunctionType,Domain,DomInfo,'below','flux');
    Qbc = [ Qbc ; QbcB ];
    boundary_points_q = [ boundary_points_q globalnr_1h(:,1)' ];
end
if bc(4) == 1
    PHIbcA = boundary_oneforms(FunctionType,Domain,DomInfo,'above','potential');
    PHIbc_ = [ PHIbc_ ; PHIbcA ];
    boundary_points_u = [ boundary_points_u globalnr_1h(:,end)' ];
else
    QbcA = boundary_oneforms(FunctionType,Domain,DomInfo,'above','flux');
    Qbc = [ Qbc ; QbcA ];
    boundary_points_q = [ boundary_points_q globalnr_1h(:,end)' ];
end

PHIbc = zeros(nr_1,1);
PHIbc(boundary_points_u) = PHIbc_;

interior_points_q = sort([ reshape(globalnr_1v((2-bc(1)):end-(1-bc(2)),:),1,[]) ...
                         reshape(globalnr_1h(:,(2-bc(3)):end-(1-bc(4))),1,[]) ]);

