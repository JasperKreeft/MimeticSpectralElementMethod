function Wbc = boundaryIntegral_assembly()

global N numElements
global nr_0 globalnr_0

disp('assembly boundary integral')

Wbce = reshape(boundaryIntegral([0 1]),(N+1)^4,1);

wbc = zeros(numElements*(N+1)^4,1);
spind0 = wbc;

for i=1:numElements

ind0 = globalnr_0(:,i);
ind_0 = (1:(N+1)^4) + (i-1)*(N+1)^4;
spind0(ind_0) = kron(ind0,ones((N+1)^2,1));
wbc(ind_0) = Wbce;

end

Wbc = sparse(spind0,spind0,wbc,nr_0,nr_0);