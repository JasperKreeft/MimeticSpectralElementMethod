function ME=MortarElement()

global Nadaptive

Na1 = Nadaptive(1);
Na2 = Nadaptive(2);

Nmin = min(Na1,Na2);

%% Mortar Left

[xi1,w1] = GLLnodes(Na1);
hm = MimeticpolyVal(xi1,Nmin,1);

W1 = diag(w1);

WM = zeros(Na1+1,Nmin+1);
for j=1:Nmin+1
    WM(:,j) = w1.*hm(j,:);
end

Z1 = W1\WM;

%% Mortar Right

[xi2,w2] = GLLnodes(Na2);
hm = MimeticpolyVal(xi2,Nmin,1);

W2 = diag(w2);

WM = zeros(Na2+1,Nmin+1);
for j=1:Nmin+1
    WM(:,j) = w2.*hm(j,:);
end

Z2 = W2\WM;

%%

ind_bc = sort( Na1+1:Na1+1:(Na1+1)^2 );
ind_in = 1:(Na1+1)^2; ind_in(ind_bc) = [];
ME1i = zeros((Na1+1)^2,Na1*(Na1+1));
ME1i(ind_in,1:Na1*(Na1+1)) = eye(Na1*(Na1+1));
ME1b = zeros((Na1+1)^2,Nmin+1);
ME1b(ind_bc,:) = Z1;


ind_bc = sort( 1:Na2+1:(Na2+1)^2 );
ind_in = 1:(Na2+1)^2; ind_in(ind_bc) = [];
ME2i = zeros((Na2+1)^2,Na2*(Na2+1));
ME2i(ind_in,1:Na2*(Na2+1)) = eye(Na2*(Na2+1));
ME2b = zeros((Na2+1)^2,Nmin+1);
ME2b(ind_bc,:) = Z2;


ME = [ ME1i zeros((Na1+1)^2,Na2*(Na2+1)) ME1b
       zeros((Na2+1)^2,Na1*(Na1+1)) ME2i ME2b ];

ME = sparse(ME);