function Wbc = boundaryIntegral_assembly()

global nr_1 globalnr_1v globalnr_1h

disp('assembly boundary integral')

Wbce = boundaryIntegral([1 2]);
        
        m1 = zeros(numElements*4*N^2*(N+1)^2,1);
        spind1r = m1;
        spind1c = m1;

        for i=1:numElements

            ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
            ind_1 = (1:4*N^2*(N+1)^2) + (i-1)*4*N^2*(N+1)^2;

spind1r(ind_1) = kron(ones(2*N*(N+1),1),ind1);
spind1c(ind_1) = kron(ind1,ones(2*N*(N+1),1));
m1(ind_1) = reshape(innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3))),4*N^2*(N+1)^2,1);

        end
        
zero = (m1==0); spind1r(zero) = []; spind1c(zero) = []; m1(zero) = [];

M  = sparse(spind1r,spind1c,m1,nr_1,nr_1);