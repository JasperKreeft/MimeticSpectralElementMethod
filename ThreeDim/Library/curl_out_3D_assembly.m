function C = curl_out_3D_assembly()

global N nr_1 nr_2
global numElements
global globalnr_1x globalnr_1y globalnr_1z
global globalnr_2x globalnr_2y globalnr_2z

disp('assembly curl')

Ce = curl_out_3D(N);

[sr,sc,ce]=find(Ce);

c = zeros(1,4*numElements*3*N^2*(N+1));
spind12_c = c;
spind12_r = c;

for i=1:numElements

    ind_12 = (1:12*N^2*(N+1)) + (i-1)*12*N^2*(N+1);

    globalnr_1 = [ globalnr_1z ; globalnr_1y ; globalnr_1x ];
    spind12_c(ind_12) = globalnr_1(sc,i);
    globalnr_2 = [ globalnr_2x ; globalnr_2y ; globalnr_2z ];
    spind12_r(ind_12) = globalnr_2(sr,i);
    c(ind_12) = ce;

end

C = sparse(spind12_r,spind12_c,c,nr_2,nr_1);

C = sign(C);