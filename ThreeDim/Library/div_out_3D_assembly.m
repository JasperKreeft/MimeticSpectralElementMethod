function D = div_out_3D_assembly()

global N nr_2 nr_3
global numElements
global globalnr_2x globalnr_2y globalnr_2z
global globalnr_3

disp('assembly div')

De = div_out_3D(N);

[sr,sc,de]=find(De);

d = zeros(1,6*numElements*N^3);
spind23_c = d;
spind23_r = d;

for i=1:numElements

    ind_23 = (1:6*N^3) + (i-1)*6*N^3;

    globalnr_2 = [ globalnr_2x ; globalnr_2y ; globalnr_2z ];
    spind23_c(ind_23) = globalnr_2(sc,i);
    spind23_r(ind_23) = globalnr_3(sr,i);
    
    d(ind_23) = de;

end

D = sparse(spind23_r,spind23_c,d,nr_3,nr_2);

D = sign(D);