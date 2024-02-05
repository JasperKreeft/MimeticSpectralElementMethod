function G = grad_out_3D_assembly()

global N nr_0 nr_1
global numElements
global globalnr_0 globalnr_1x globalnr_1y globalnr_1z

disp('assembly gradient')

Ge = grad_out_3D(N);

[sr,sc,ge]=find(Ge);

g = zeros(1,2*numElements*3*N*(N+1)^2);
spind01_c = g;
spind01_r = g;

for i=1:numElements

    ind_01 = (1:6*N*(N+1)^2) + (i-1)*6*N*(N+1)^2;

    spind01_c(ind_01) = globalnr_0(sc,i);
    globalnr_1 = [ globalnr_1z ; globalnr_1y ; globalnr_1x ];
    spind01_r(ind_01) = globalnr_1(sr,i);
    g(ind_01) = ge;

end

G = sparse(spind01_r,spind01_c,g,nr_1,nr_0);

G = sign(G);