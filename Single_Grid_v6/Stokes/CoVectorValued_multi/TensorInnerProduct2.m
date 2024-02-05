function T2 = TensorInnerProduct2()
% Tensor inner product for volume forms

global N xi w

[xiw,ww] = Gnodes(N);

[~,e] = MimeticpolyVal(xi,N,1);
[~,ew] = MimeticpolyVal(xiw,N,3);

Ix_pik = zeros(N+1);
for i=1:N+1
    for k=1:N+1
        Ix_pik(k,i) = sum(ww.*ew(i,:).*ew(k,:));
    end
end

Ix_qjl = zeros(N);
for j=1:N
    for l=1:N
        Ix_qjl(l,j) = sum(w.*e(j,:).*e(l,:));
    end
end

Ix_x = kron(ones(N),Ix_pik);

Ix_y = kron(Ix_qjl,ones(N+1));

Ix = Ix_x.*Ix_y;

% Iy = Ix;

T2 = kron(speye(2),Ix);