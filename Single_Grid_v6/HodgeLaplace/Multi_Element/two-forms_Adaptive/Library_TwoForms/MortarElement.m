function ME=MortarElement()

global Nadaptive Nmax

% Na1<=Na2

Na1 = Nadaptive(1);
Na2 = Nadaptive(2);

Nmin = min(Na1,Na2);

%% Mortar Left

[xi1,w1] = GLLnodes(Na1);

[~,e1] = MimeticpolyVal(xi1,Na1,1);
[~,em] = MimeticpolyVal(xi1,Nmin,1);


W1 = zeros(Na1);
for l=1:Na1
    for j=1:Na1
        W1(l,j) = sum(w1.*e1(j,:).*e1(l,:));
    end
end

WM = zeros(Na1,Nmin);
for l=1:Na1
    for k=1:Nmin
        WM(l,k) = sum(w1.*em(k,:).*e1(l,:));
    end
end

Z1 = W1\WM;

%% Mortar Right

[xi2,w2] = GLLnodes(Na2);

[~,e2] = MimeticpolyVal(xi2,Na2,1);
[~,em] = MimeticpolyVal(xi2,Nmin,1);


W2 = zeros(Na2);
for l=1:Na2
    for j=1:Na2
        W2(l,j) = sum(w2.*e2(j,:).*e2(l,:));
    end
end

WM = zeros(Na2,Nmin);
for l=1:Na2
    for k=1:Nmin
        WM(l,k) = sum(w2.*em(k,:).*e2(l,:));
    end
end

Z2 = W2\WM;

%%


indv = Na1*(Na1+1)+Na2*(Na2+1);
indN = Na1^2+Na1*(Na1+1);
indM = Na2^2+Na2*(Na2+1);

ME1 = zeros(2*Na1*(Na1+1),indN+indM+Nmin);
ind_bc = Na1+1:Na1+1:Na1*(Na1+1);
ind_in = 1:Na1*(Na1+1); ind_in(ind_bc) = [];
ME1(ind_in,1:Na1^2) = eye(Na1^2);
ME1(ind_bc,indN+indM+(1:Nmin)) = Z1;
ME1(Na1*(Na1+1)+(1:Na1*(Na1+1)),Na1^2+(1:Na1*(Na1+1))) = eye(Na1*(Na1+1));

ME2 = zeros(2*Na2*(Na2+1),indN+indM+Nmin);
ind_bc = 1:Na2+1:Na2*(Na2+1);
ind_in = 1:Na2*(Na2+1); ind_in(ind_bc) = [];
ME2(ind_in,indN+(1:Na2^2)) = eye(Na2^2);
ME2(ind_bc,indN+indM+(1:Nmin)) = Z2;
ME2(Na2*(Na2+1)+(1:Na2*(Na2+1)),indN+Na2^2+(1:Na2*(Na2+1))) = eye(Na2*(Na2+1));

MEflux = [ ME1
           ME2 ];

ME = [ MEflux zeros(2*Na1*(Na1+1)+2*Na2*(Na2+1),Na1^2+Na2^2)
       zeros(Na1^2+Na2^2,indN+indM+Nmin) eye(Na1^2+Na2^2) ];

ME = sparse(ME);