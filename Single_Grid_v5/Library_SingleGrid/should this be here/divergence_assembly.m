function D = divergence_assembly()

global N nr_1 nr_2
global numElements
global globalnr_2 globalnr_1v globalnr_1h

disp('assembly divergence')

De = div(N);

[sr,sc,de]=find(De);

d = zeros(1,4*numElements*N^2);
spind12_c = d;
spind12_r = d;

for i=1:numElements

    ind_12 = (1:4*N^2) + (i-1)*4*N^2;
    
    globalnr_1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
    spind12_c(ind_12) = globalnr_1(sc);
    
    spind12_r(ind_12) = globalnr_2(sr,i);
    d(ind_12) = de;

end

D = sparse(spind12_r,spind12_c,d,nr_2,nr_1);

D = sign(D);


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