function M = innerproduct_assembly(form,Mesh)

global N
global numElements

switch form

    case 0 % zero-forms

        disp('assembly innerproduct zero-forms')
        global nr_0 globalnr_0 
        Ne = (N+1)^2; nr = nr_0;
        globalnr = globalnr_0;

    case 1 % one-forms

        disp('assembly innerproduct one-forms')
        global nr_1 globalnr_1v globalnr_1h
        Ne = 2*N*(N+1); nr = nr_1;
        globalnr = [ globalnr_1v ; globalnr_1h ];

    case 2 % two-forms

        disp('assembly innerproduct two-forms')
        global nr_2 globalnr_2
        Ne = N^2; nr = nr_2;
        globalnr = globalnr_2;

end

switch form

    case 0 % zero-forms

        m = zeros(numElements*Ne,1); % if form==0; Ne = sqrt(Ne); end % WHY ??????
        spind_r = m;
        spind_c = m;

        for i=1:numElements

            ind = globalnr(:,i);
            ind_ = (1:Ne) + (i-1)*Ne;

            spind(ind_) = ind;
            m(ind_) = diag(innerproduct(0,Mesh.J(:,i)));

        end


    otherwise % one- or two-forms

        m = zeros(numElements*Ne^2,1);
        spind_r = m;
        spind_c = m;

        for i=1:numElements

            ind = globalnr(:,i);

            ind_ = (1:Ne^2) + (i-1)*Ne^2;

            spind_r(ind_) = kron(ones(Ne,1),ind);
            spind_c(ind_) = kron(ind,ones(Ne,1));

            if form==1
                me = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
            elseif form==2
                me = innerproduct(2,Mesh.J(:,i));
            end

            m(ind_) = reshape(me,Ne^2,1);

        end

end

if form==1
    zero = (m==0); spind_r(zero) = []; spind_c(zero) = []; m(zero) = [];
end

M  = sparse(spind_r,spind_c,m,nr,nr);


%#ok<*TLEV>