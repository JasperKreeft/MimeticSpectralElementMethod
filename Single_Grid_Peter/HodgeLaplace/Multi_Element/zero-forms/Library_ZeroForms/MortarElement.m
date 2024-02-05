function ME=MortarElement()

global Nadaptive Nmax

% Na1<=Na2

Na1 = Nadaptive(1);
Na2 = Nadaptive(2);

[xi,w] = GLLnodes(Nmax);

h1 = MimeticpolyVal(xi,Na1,1);

W1 = diag(w);

W2 = zeros(Nmax+1,Na1+1);
for q=1:Nmax+1
    for j=1:Na1+1
        W2(q,j) = w(q)*h1(j,q);
    end
end

Z = W1\W2;

ME1 = zeros((Na1+1)^2,Na1*(Na1+1)+Na2*(Na2+1)+Na1+1);
ind =[]; for i=1:Na1+1; ind = [ ind (i-1)*(Na1+1)+(1:Na1) ]; end
ME1(ind,1:Na1*(Na1+1)) = eye(Na1*(Na1+1));
ME1(Na1+1:Na1+1:(Na1+1)^2,Na1*(Na1+1)+Na2*(Na2+1)+(1:Na1+1)) = eye(Na1+1);


ind_bc = sort( 1:Na2+1:(Na2+1)^2 );

ind_in = 1:(Na2+1)^2; ind_in(ind_bc) = [];


ME2 = zeros((Na2+1)^2,Na2*(Na2+1)+(Na1+1));

ME2(ind_in,1:length(ind_in)) = eye(length(ind_in));
ME2(ind_bc,Na2*(Na2+1)+(1:Na1+1)) = Z;

ME = [ ME1
       zeros((Na2+1)^2,Na1*(Na1+1)) ME2 ];




