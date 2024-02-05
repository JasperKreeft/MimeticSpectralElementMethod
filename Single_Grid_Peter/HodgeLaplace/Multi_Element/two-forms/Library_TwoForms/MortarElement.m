function ME=MortarElement()

global Nadaptive Nmax

% Na1<=Na2

Na1 = Nadaptive(1);
Na2 = Nadaptive(2);

[xi,w] = GLLnodes(Nmax);

[~,e1] = MimeticpolyVal(xi,Na1,1);
[~,e2] = MimeticpolyVal(xi,Na2,1);


W1 = zeros(Nmax,Na2);
for l=1:Nmax
    for j=1:Na2
        W1(l,j) = sum(w.*e2(j,:).*e2(l,:));
    end
end

W2 = zeros(Nmax,Na1);
for l=1:Nmax
    for k=1:Na1
        W2(l,k) = sum(w.*e1(k,:).*e2(l,:));
    end
end

Z = W1\W2;


indv = Na1*(Na1+1)+Na2*(Na2+1);
indN = Na1^2+Na1*(Na1+1);
indM = Na2^2+Na2*(Na2+1);

ME1 = zeros(2*Na1*(Na1+1),Na1^2+Na2^2+Na1+indv);
ind =[]; for i=1:Na1; ind = [ ind (i-1)*(Na1+1)+(1:Na1) ]; end
ME1(ind,1:Na1^2) = eye(Na1^2);
ME1(Na1+1:Na1+1:Na1*(Na1+1),Na1^2+Na2^2+indv+(1:Na1)) = eye(Na1);
ME1(Na1*(Na1+1)+(1:Na1*(Na1+1)),Na1^2+(1:Na1*(Na1+1))) = eye(Na1*(Na1+1));


ind_bc = sort( 1:Na2+1:Na2*(Na2+1) );

ind_in = 1:Na2*(Na2+1); ind_in(ind_bc) = [];

ME2 = zeros(8,Na1^2+Na2^2+Na1+indv);

ME2(ind_in,indN+(1:Na2^2)) = eye(Na2^2);
ME2(ind_bc,indN+Na2^2+Na2*(Na2+1)+(1:Na1)) = Z;
ME2(Na2*(Na2+1)+(1:Na2*(Na2+1)),indN+Na2^2+(1:Na2*(Na2+1))) = eye(Na2*(Na2+1));

MEflux = [ ME1
           ME2 ];

ME = [ MEflux zeros(2*Na1*(Na1+1)+2*Na2*(Na2+1),Na1^2+Na2^2)
       zeros(Na1^2+Na2^2,indN+indM+Na1) eye(Na1^2+Na2^2) ];

ME = sparse(ME);