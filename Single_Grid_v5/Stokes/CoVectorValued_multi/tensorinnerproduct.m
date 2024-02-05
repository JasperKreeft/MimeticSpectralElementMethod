function T = tensorinnerproduct()

global N w e

% [xi_,w_] = GLLnodes(N+2);
[xi_,w_] = Gnodes(N+2);
[hw,~] = MimeticpolyVal(xi_,N,3);

XXYY = zeros(N*(N+2));
% xx - yy
for k=1:N+2
    for l=1:N
        kl = k+(l-1)*(N+2);
        for i=1:N+2
            for j=1:N
                ij = i+(j-1)*(N+2);
                
                XXYY(kl,ij) = sum(w_.*hw(i,:).*hw(k,:))*sum(w.*e(j,:).*e(l,:));
                
            end
        end
    end
end


% [xi_,w_] = GLLnodes(N+1);
[xi_,w_] = Gnodes(N+1);
[~,ew] = MimeticpolyVal(xi_,N,3);

% yx - xy
YXXY = zeros((N+1)^2);
for k=1:N+1
    for q=1:N+1
        kq = k+(q-1)*(N+1);
        for i=1:N+1
            iq = i+(q-1)*(N+1);

            YXXY(kq,iq) = w(q)*sum(w_.*ew(i,:).*ew(k,:));

        end
    end
end


[rows1,cols1,vals1] = find(XXYY);
[rows2,cols2,vals2] = find(YXXY);

rows = [ rows1 ; rows2+N*(N+2) ];
rows = [ rows ; rows + N*(N+2)+(N+1)^2 ];
cols = [ cols1 ; cols2+N*(N+2) ];
cols = [ cols ; cols + N*(N+2)+(N+1)^2 ];
vals = [ vals1; vals2 ; vals1; vals2 ];

T = sparse(rows,cols,vals,2*(N*(N+2)+(N+1)^2),2*(N*(N+2)+(N+1)^2));
