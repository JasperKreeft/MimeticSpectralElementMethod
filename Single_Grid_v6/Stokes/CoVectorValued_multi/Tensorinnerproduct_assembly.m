function T = Tensorinnerproduct_assembly()

% ALLEEN GELDIG VOOR MEERDERE STANDAARD ELEMENTEN

global N
global numElements
global nr_11
global globalnr_11_xx globalnr_11_yx globalnr_11_yy globalnr_11_xy

disp('assembly innerproduct covector-valued (n-1)-forms')

nr_11_e   = 2*N*(N+2)+2*(N+1)^2;
% nr_coeff_e = 2*N^2*(N+2)^2+2*(N+1)^4;

t1 = zeros(numElements*nr_11_e^2,1);
spind1r = t1;
spind1c = t1;

Te = tensorinnerproduct();

te = zeros(nr_11_e^2,1);
[sr,sc,t_] = find(Te);
te(sr+(sc-1)*nr_11_e) = t_;

for i=1:numElements

    ind1  = [ globalnr_11_xx(:,i) ; globalnr_11_yx(:,i) ; globalnr_11_yy(:,i) ; globalnr_11_xy(:,i) ];
    ind_1 = (1:nr_11_e^2) + (i-1)*nr_11_e^2;

    spind1r(ind_1) = kron(ones(nr_11_e,1),ind1);
    spind1c(ind_1) = kron(ind1,ones(nr_11_e,1));
    t1(ind_1) = te; % t1(ind_1) + 

end

zero = (t1==0); spind1r(zero) = []; spind1c(zero) = []; t1(zero) = [];

T  = sparse(spind1r,spind1c,t1,nr_11,nr_11);
