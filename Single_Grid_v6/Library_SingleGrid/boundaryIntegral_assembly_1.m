function Wbc = boundaryIntegral_assembly_1()

global N numElements
global nr_1 globalnr_1v globalnr_1h

disp('assembly boundary integral')

Wbce = reshape(boundaryIntegral([1 2]),4*N^2*(N+1)^2,1);

wbc = zeros(numElements*4*N^2*(N+1)^2,1);
spind1r = wbc;
spind1c = wbc;

for i=1:numElements

ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
ind_1 = (1:4*N^2*(N+1)^2) + (i-1)*4*N^2*(N+1)^2;
spind1r(ind_1) = kron(ones(2*N*(N+1),1),ind1);
spind1c(ind_1) = kron(ind1,ones(2*N*(N+1),1));
wbc(ind_1) = Wbce;

end

Wbc = sparse(spind1r,spind1c,wbc,nr_1,nr_1);



%% Full-memory
% 
% Wbc1 = zeros(nr_1);
% for i=1:numElements
% ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
% Wbc1(ind1,ind1) = Wbc1(ind1,ind1) + boundaryIntegral([1 2]);
% end
% Wbc1 = sparse(Wbc1);