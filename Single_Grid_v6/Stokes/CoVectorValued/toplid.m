function [Wbc,UVbc]=toplid()

global N

P=N+1;
[xip,wp] = Gnodes(P);
[hegp,eegp] = MimeticpolyVal(xip,N,3);

A = zeros((N+2)*(N+1),1); row = A; col = A;
for i=1:N+2
    for k=1:N+1
        ki = k+(i-1)*(N+1);
        A(ki,1) = sum(wp.*hegp(i,:).*eegp(k,:));
%         row(ki,1) = (N+2)*N + N*(N+1)+k;
        row(ki,1) = (N+2)*N + k*(N+1);
        col(ki,1) = (N+2)*(N+1) + i*(N+1);
    end
end

Wbc = sparse(row,col,A,2*(N*(N+2)+(N+1)^2),2*(N+2)*(N+1));



%%%

v = ones(N+2,1);

UVbc = zeros(2*(N+2)*(N+1),1);

UVbc((N+2)*(N+1)+(1:N+2)*(N+1),1) = v;