function ME=MortarElement()

global Nadaptive numRows numColumns numElements
global Nmax nr_0

% Na1 = Nadaptive(1);
% Na2 = Nadaptive(2);

Na = zeros(numColumns*(numRows-1)+numRows*(numColumns-1),2);
for c=1:numColumns-1
    rc = (1:numRows)+(c-1)*numRows;
    Na(rc,:) = Nadaptive(:,c:c+1);
end
for r=1:numRows-1
    rc = (1:numColumns) + (r-1)*numColumns + numRows*(numColumns-1);
    Na(rc,:) = Nadaptive(r:r+1,:)';
end
Na = unique(Na,'rows');


Zlocal = cell(Nmax,Nmax);
for i=1:size(Na,1)
    
Namin = min(Na(i,:));
Na1 = Na(i,1);
Na2 = Na(i,2);

%% Mortar Left

[xi1,w1] = GLLnodes(Na1);
hm = MimeticpolyVal(xi1,Namin,1);

W1 = diag(w1);

WM = zeros(Na1+1,Namin+1);
for j=1:Namin+1
    WM(:,j) = w1.*hm(j,:);
end

Z1 = W1\WM;

%% Mortar Right

[xi2,w2] = GLLnodes(Na2);
hm = MimeticpolyVal(xi2,Namin,1);

W2 = diag(w2);

WM = zeros(Na2+1,Namin+1);
for j=1:Namin+1
    WM(:,j) = w2.*hm(j,:);
end

Z2 = W2\WM;

%%

ind_bc = Na1+1:Na1+1:(Na1+1)^2;
ind_in = 1:(Na1+1)^2; ind_in(ind_bc) = [];
ME1i = zeros((Na1+1)^2,Na1*(Na1+1));
ME1i(ind_in,1:Na1*(Na1+1)) = eye(Na1*(Na1+1));
ME1b = zeros((Na1+1)^2,Namin+1);
ME1b(ind_bc,:) = Z1;


ind_bc = 1:Na2+1:(Na2+1)^2;
ind_in = 1:(Na2+1)^2; ind_in(ind_bc) = [];
ME2i = zeros((Na2+1)^2,Na2*(Na2+1));
ME2i(ind_in,1:Na2*(Na2+1)) = eye(Na2*(Na2+1));
ME2b = zeros((Na2+1)^2,Namin+1);
ME2b(ind_bc,:) = Z2;


ME = [ ME1i zeros((Na1+1)^2,Na2*(Na2+1)) ME1b
       zeros((Na2+1)^2,Na1*(Na1+1)) ME2i ME2b ];

Zlocal{Na1,Na2} = sparse(ME);
keyboard
end

%%%%%%%%%%%%%%%%%%%%

[cnectN cnectL] = elConnect(numColumns,numRows);

nr_0_in_element = (Nadaptive+1).^2;
for i=1:numElements
    neigb = cnectN(i,:); neigb(neigb==0) = [];
    for j=neigb
        dif_N = Nadaptive(i) - Nadaptive(j);
        nr_0_in_element(i) = nr_0_in_element(i) - (dif_N>0)*dif_N;
    end
end

ind = nr_0 - sum(max(Na'));
Z = zeros(nr_0,ind);
% for i=1:numElements
%     ind1 = globalnr_0(1:(Nadaptive(i)^2+1),i);
%     ind2 = (1:nr_0_in_element(i)) + nr_0_in_element(1:i-1);
%     
%     Z(ind1,ind2) = Zlocal{}
% 
% 
% 




