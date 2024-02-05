function M = innerproduct_3D_assembly(form,Mesh)

global N
global numElements

tic
switch form


    case 0 % zero-forms

        global nr_0 globalnr_0 

        disp('assembly innerproduct zero-forms')

        m0 = zeros(numElements*(N+1)^6,1);
        spind0 = m0;

        for i=1:numElements

            ind0 = globalnr_0(:,i);
            ind_0 = (1:(N+1)^6) + (i-1)*(N+1)^6;

            spind0(ind_0) = kron(ind0,ones((N+1)^3,1));
%             m0(ind_0) = reshape(innerproduct(0,Mesh.J(:,i)),(N+1)^6,1);
            m0(ind_0) = reshape(innerproduct_3D(0),(N+1)^6,1);

        end

        M  = sparse(spind0,spind0,m0,nr_0,nr_0);


    case 1 % one-forms

        global nr_1 globalnr_1x globalnr_1y globalnr_1z

        disp('assembly innerproduct one-forms')

        m1 = zeros(numElements*9*N^2*(N+1)^4,1);
        spind1r = m1;
        spind1c = m1;

        for i=1:numElements

            ind1 = [ globalnr_1z(:,i) ; globalnr_1y(:,i) ; globalnr_1x(:,i) ];
            ind_1 = (1:9*N^2*(N+1)^4) + (i-1)*9*N^2*(N+1)^4;

spind1r(ind_1) = kron(ones(3*N*(N+1)^2,1),ind1);
spind1c(ind_1) = kron(ind1,ones(3*N*(N+1)^2,1));
m1(ind_1) = reshape(innerproduct_3D(1),9*N^2*(N+1)^4,1);

        end

        nz = (m1==0); spind1r(nz) = []; spind1c(nz) = []; m1(nz) = [];
        M  = sparse(spind1r,spind1c,m1,nr_1,nr_1);


    case 2 % one-forms

        global nr_2 globalnr_2x globalnr_2y globalnr_2z

        disp('assembly innerproduct two-forms')

        m2 = zeros(numElements*9*N^4*(N+1)^2,1);
        spind2r = m2;
        spind2c = m2;

        for i=1:numElements

            ind2 = [ globalnr_2x(:,i) ; globalnr_2y(:,i) ; globalnr_2z(:,i) ];
            ind_2 = (1:9*N^4*(N+1)^2) + (i-1)*9*N^4*(N+1)^2;

spind2r(ind_2) = kron(ones(3*N^2*(N+1),1),ind2);
spind2c(ind_2) = kron(ind2,ones(3*N^2*(N+1),1));
m2(ind_2) = reshape(innerproduct_3D(2),9*N^4*(N+1)^2,1);

        end

        nz = (m2==0); spind2r(nz) = []; spind2c(nz) = []; m2(nz) = [];
        M  = sparse(spind2r,spind2c,m2,nr_2,nr_2);



    case 3 % two-forms

        global nr_3 globalnr_3

        disp('assembly innerproduct three-forms')

        m3 = zeros(numElements*N^6,1);
        spind3r = m3;
        spind3c = m3;

        for i=1:numElements

            ind3 = globalnr_3(:,i);
            ind_3 = (1:N^6) + (i-1)*N^6;

            spind3r(ind_3) = kron(ones(N^3,1),ind3);
            spind3c(ind_3) = kron(ind3,ones(N^3,1));
            m3(ind_3) = reshape(innerproduct_3D(3),N^6,1);

        end

        M  = sparse(spind3r,spind3c,m3,nr_3,nr_3);

end
toc

%#ok<*TLEV>