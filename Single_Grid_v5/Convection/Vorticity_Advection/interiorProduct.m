function vE01 = interiorProduct(U,V,Exe,Eye,he,eh,Mesh)

disp('Assemble interior product')

global N numElements
global globalnr_0 globalnr_1v globalnr_1h
global nr_0 nr_1

vE01 = zeros(nr_0,nr_1);       
for i=1:numElements

ub = (he*U(:,i))./Mesh.J(:,i);
vb = (eh*V(:,i))./Mesh.J(:,i);

Axc = zeros((N+1)^2,N+1);
Ayc = zeros((N+1)^2,N+1);
for j=1:N+1
    ind = (j-1)*(N+1)+(1:N+1);
    Axc(ind,j) = ub(ind);
    Ayc(ind,:) = diag(vb(ind));
end
Axe = sparse(repmat(Axc,1,N));
Aye = sparse(repmat(Ayc,1,N));

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];

vE01(ind0,ind1) = vE01(ind0,ind1) + [ -Aye.*Eye Axe.*Exe ];

end

ind_element_boundary = [1:N+1 N+2:N+1:N*(N+1) 2*(N+1):N+1:N*(N+1) N*(N+1)+1:(N+1)^2];
ind_all_el_bcs = unique(reshape(globalnr_0(ind_element_boundary,:),[],1));
vE01(ind_all_el_bcs,:) = vE01(ind_all_el_bcs,:)/2; % central discretization

indCornerPoints = unique(globalnr_0([1 N+1 N*(N+1)+1 (N+1)^2],:));
vE01(indCornerPoints,:) = vE01(indCornerPoints,:)/2; % central discretization

vE01 = sparse(vE01);