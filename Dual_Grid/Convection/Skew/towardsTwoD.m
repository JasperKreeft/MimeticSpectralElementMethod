clear all
close all
clc

global N w

N = 2;

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);


Ex = kron(eye(N),e');
Ey = zeros(N*(N+1),N*N);
for i=1:N
    Ey(1:N*(N+1),(i-1)*N+(1:N)) = kron(eye(N),e(i,:)');
end
E = [ Ex; Ey];
        
D = div(N);

M2 = innerproduct_twoforms(e);

K1 = M2*D*E;
K2 = M2;



t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%%

    phinew = timemarching(K2,K1,dt,phi,'ESDIRK3');
    
    pphinew = phinew'*ee;

    phidx = phinew./diff(xi)';
    plot([xi(2:N+1) ; xi(1:N)],[phidx phidx]','g')
    plot(xx,pphinew)
%     ylim([-.2 1])
    pause(0.05)
    hold off
    phi = phinew;
    
end