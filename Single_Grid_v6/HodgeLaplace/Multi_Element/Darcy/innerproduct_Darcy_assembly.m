function M = innerproduct_Darcy_assembly(form,Mesh)

global N
global numElements

switch form



    case 0 % zero-forms

        global nr_0 globalnr_0 

        disp('assembly innerproduct zero-forms')

        m0 = zeros(numElements*(N+1)^2,1);
        spind0 = m0;

        for i=1:numElements

            ind0 = globalnr_0(:,i);
            ind_0 = (1:(N+1)^2) + (i-1)*(N+1)^2;

            spind0(ind_0) = ind0;
            m0(ind_0) = diag(innerproduct(0,Mesh.J(:,i)));

        end

        M  = sparse(spind0,spind0,m0,nr_0,nr_0);


    case 1 % one-forms

        global nr_1 globalnr_1v globalnr_1h
        global i

        disp('assembly innerproduct one-forms')

        m1 = zeros(numElements*4*N^2*(N+1)^2,1);
        spind1r = m1;
        spind1c = m1;

        for i=1:numElements

            ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
            ind_1 = (1:4*N^2*(N+1)^2) + (i-1)*4*N^2*(N+1)^2;

spind1r(ind_1) = kron(ones(2*N*(N+1),1),ind1);
spind1c(ind_1) = kron(ind1,ones(2*N*(N+1),1));

m1e = innerproduct_Darcy(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)),Mesh.X(:,i),Mesh.Y(:,i));
m1(ind_1) = reshape(m1e,4*N^2*(N+1)^2,1);

        end
        
zero = (m1==0); spind1r(zero) = []; spind1c(zero) = []; m1(zero) = [];

M  = sparse(spind1r,spind1c,m1,nr_1,nr_1);



    case 2 % two-forms

        global nr_2 globalnr_2

        disp('assembly innerproduct two-forms')

        m2 = zeros(numElements*N^4,1);
        spind2r = m2;
        spind2c = m2;

        for i=1:numElements

            ind2 = globalnr_2(:,i);
            ind_2 = (1:N^4) + (i-1)*N^4;

            spind2r(ind_2) = kron(ones(N^2,1),ind2);
            spind2c(ind_2) = kron(ind2,ones(N^2,1));
            m2(ind_2) = reshape(innerproduct(2,Mesh.J(:,i)),N^4,1);

        end

        M  = sparse(spind2r,spind2c,m2,nr_2,nr_2);

end


%#ok<*TLEV>