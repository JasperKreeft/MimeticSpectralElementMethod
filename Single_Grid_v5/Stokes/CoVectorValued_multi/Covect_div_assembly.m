function CD = Covect_div_assembly()

global N nr_11 nr_21
global numElements
global globalnr_11_xx globalnr_11_yx globalnr_11_yy globalnr_11_xy
global globalnr_21_x globalnr_21_y

disp('assembly covector divergence')

CDe = covect_div(N);

[sr,sc,cde]=find(CDe);

cd = zeros(1,4*numElements*2*N*(N+1));
spind12_c = cd;
spind12_r = cd;

for i=1:numElements

    ind_12 = (1:4*2*N*(N+1)) + (i-1)*4*2*N*(N+1);
    
    globalnr_11 = [ globalnr_11_xx(:,i) ; globalnr_11_yx(:,i) ; globalnr_11_yy(:,i) ; globalnr_11_xy(:,i) ];
    spind12_c(ind_12) = globalnr_11(sc);
    
    globalnr_21 = [ globalnr_21_x(:,i) ; globalnr_21_y(:,i) ];
    spind12_r(ind_12) = globalnr_21(sr);
    cd(ind_12) = cde;

end

CD = sparse(spind12_r,spind12_c,cd,nr_21,nr_11);

CD = sign(CD);


%% full-memory version
% D  = zeros(nr_2,nr_1);
% 
% for i=1:numElements
% 
% ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
% ind2 = globalnr_2(:,i);
% 
% D(ind2,ind1) = De;
% 
% end
% 
% D = sparse(D);