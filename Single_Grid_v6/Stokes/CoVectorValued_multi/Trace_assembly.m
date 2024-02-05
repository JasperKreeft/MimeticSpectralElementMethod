function Tr = Trace_assembly()

% ALLEEN GELDIG VOOR MEERDERE STANDAARD ELEMENTEN

global N
global numElements
global nr_1 nr_21
global globalnr_1v globalnr_1h
global globalnr_21_x globalnr_21_y

disp('assembly Trace from covector-valued (n)-forms to real-valued (n-1)-forms')


tr = zeros(numElements*4*N^2*(N+1)^2,1);
spind211r = tr;
spind211c = tr;

Tre = Trace();

[sr,sc,tre] = find(Tre);

for i=1:numElements

    globalnr_1  = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
    globalnr_21 = [ globalnr_21_x(:,i) ; globalnr_21_y(:,i) ];
    ind_211 = (1:2*N*(N+1)^2) + (i-1)*2*N*(N+1)^2;

    spind211r(ind_211) = globalnr_1(sr);
    spind211c(ind_211) = globalnr_21(sc);
    tr(ind_211) = tre; %tr(ind_211) + 

end

zero = (tr==0); spind211r(zero) = []; spind211c(zero) = []; tr(zero) = [];

Tr  = sparse(spind211r,spind211c,tr,nr_1,nr_21);
