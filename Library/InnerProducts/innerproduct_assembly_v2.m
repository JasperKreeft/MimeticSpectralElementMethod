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

m = zeros(numElements*Ne^2,1);
spind_r = m;
spind_c = m;

    for i=1:numElements

%         if form==1
%             m(ind_) = reshape(innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3))),Ne^2,1);
%         else
%             m(ind_) = reshape(innerproduct(form,Mesh.J(:,i)),Ne^2,1);
%         end
%         ind = globalnr(:,i);
%         ind_ = (1:Ne^2) + (i-1)*Ne^2;
% 
%         spind_r(ind_) = kron(ones(Ne,1),ind);
%         spind_c(ind_) = kron(ind,ones(Ne,1));

        if form==1
            me = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
        else
            me = innerproduct(form,Mesh.J(:,i));
        end
        [sr,sc,m_e]=find(me);

%         nz = (sr==0); sr(nz) = []; sc(nz) = []; m_e(nz) = [];
        
        ind = globalnr(:,i);
        ss = size(sr,1);
        ind_ = (1:ss) + (i-1)*ss;

        spind_r(ind_) = ind(sr);
        spind_c(ind_) = ind(sc);
        
        m(ind_) = m_e;

    end
    
    nz = (spind_r==0); spind_r(nz) = []; spind_c(nz) = []; m(nz) = [];

    M  = sparse(spind_r,spind_c,m,nr,nr);

end


%#ok<*TLEV>