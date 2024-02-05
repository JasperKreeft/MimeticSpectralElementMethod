
addpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
addpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Library_HodgeLaplace')
addpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Multi_Element/Library_MultiElement')

global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2
global xi

xi = GLLnodes(N);

numbering('square')

Qxi  = Q(globalnr_1v);
Qeta = Q(globalnr_1h);

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square('SinDeformGrid',cc);

[qx,qy,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);

[Qxi_interp,Qeta_interp] = fluxValue('sine','SinDeformGrid',cc);
[qx_interp,qy_interp,q_interp] = reconstruct(1,Qxi_interp,Qeta_interp,hp,ep,Meshp);

[qx_ex,qy_ex] = exact_solution(Meshp.X,Meshp.Y,'sine','one');






er = er+1;

if isempty(numElements)
    numElements = numRows*numColumns;
end

for i=1:numElements

JW = Meshp.J(:,i).*Meshp.W;

% %% Flux error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorL2_q(er)        = errorL2_q(er)+sum(sum((qx_ex(:,i)-qx(:,i)).^2.*JW));
if exist('qx_interp','var')
errorL2_q_interp(er) = errorL2_q_interp(er)+sum(sum((qx_ex(:,i)-qx_interp(:,i)).^2.*JW));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorL2_q(er)        = sqrt(abs(errorL2_q(er)));
errorL2_q_interp(er) = sqrt(abs(errorL2_q_interp(er)));



plotten


rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Library_HodgeLaplace')
rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Multi_Element/Library_MultiElement')




pause