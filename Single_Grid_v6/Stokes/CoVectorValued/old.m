clear all
close all
clc

N=2;

PQ = N+2;

[xi,w] = GLLnodes(N);
[xipq,wpq] = GLLnodes(PQ);

[h,e] = MimeticpolyVal(xi,N,1);
[hw,~] = MimeticpolyVal(xipq,N,3);

XXYY = zeros(N*(N+2));
% xx - yy
for k=1:N+2
    for l=1:N
        kl = k+(l-1)*(N+2);
        for i=1:N+2
            for j=1:N
                ij = i+(j-1)*(N+2);
                
                XXYY(kl,ij) = sum(wpq.*hw(i,:).*hw(k,:))*sum(w.*e(j,:).*e(l,:));
                
            end
        end
    end
end



PQ = N+2;

% [xipq,wpq] = GLLnodes(PQ);
[xipq,wpq] = Gnodes(PQ);
[hw,ew] = MimeticpolyVal(xipq,N,3);


% yy
for k=1:N
    for l=1:N+2
        kl = l+(k-1)*(N+2);
        for i=1:N
            for j=1:N+2
                ij = j+(i-1)*(N+2);
                
                YY(kl,ij) = sum(w.*e(i,:).*e(k,:))*sum(wpq.*hw(j,:).*hw(l,:));
                
            end
        end
    end
end

XXYY
YY
abs(XXYY-YY)



PQ = N+1;

[xipq,wpq] = GLLnodes(PQ);
% [xipq,wpq] = Gnodes(PQ);
[~,ew] = MimeticpolyVal(xipq,N,3);

% yx
for k=1:N+1
    for q=1:N+1
        kq = q+(k-1)*(N+1);
        for i=1:N+1
            iq = q+(i-1)*(N+1);
            
            YX(kq,iq) = w(q)*sum(wpq.*ew(i,:).*ew(k,:));
            
        end
    end
end


% xy
for p=1:N+1
    for l=1:N+1
        pl = p+(l-1)*(N+1);
        for j=1:N+1
            pj = p+(j-1)*(N+1);
            
            XY(pl,pj) = w(p)*sum(wpq.*ew(j,:).*ew(l,:));
            
        end
    end
end

clc
YX
XY
abs(YX-XY)