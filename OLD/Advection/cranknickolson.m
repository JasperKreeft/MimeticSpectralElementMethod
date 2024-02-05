function [X_CN] = cranknickolson(dT,Xinit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CN / trapeziodal with pseudo time stepping

NT = 17;
dt = dT/(NT-1);

dtau = dt;

t    = zeros(1,NT);
X_CN = zeros(2,NT);

X_CN(:,1) = Xinit;

for n=1:NT-1
    t(n+1) = t(n)+dt;
    [A_n,B_n] = vector_elements(X_CN(:,n));

    Xn_new = X_CN(:,n);
    Xn_old = [10;10];
    while abs(max(Xn_new-Xn_old)./Xn_old)>1e-12
        Xn_old = Xn_new;
        [A_np1,B_np1] = vector_elements(Xn_old);
        Xn_new = Xn_old-dtau*((Xn_old-X_CN(:,n))/dt-1/2*(A_np1*Xn_old+B_np1)-1/2*(A_n*X_CN(:,n)+B_n));
    end
    X_CN(:,n+1) = Xn_new;
end
plot(X_CN(2,:),X_CN(1,:),'x-r','linewidth',2)