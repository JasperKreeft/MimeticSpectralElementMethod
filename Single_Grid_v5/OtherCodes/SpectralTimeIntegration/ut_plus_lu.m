clear all
clf
clc

N = 12;

tbegin = 0;
tend   = 1;
dt = tend-tbegin;

phi0 = 1;

lambda = -10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 100000;
tt = linspace(-1,1,n);
hh = LagrangeVal(tt,N,1);

tt = (tbegin+tend)/2+(tend-tbegin)/2*tt;

phi_ex = phi0*exp(lambda*tt);

subplot(1,2,1)
plot(tt,phi_ex,'r')
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[tGLL] = GLLnodes(N);

difft = diff(tGLL)';

[sG,wG] = Gnodes(ceil((N+1)/2));

A = diag(ones(N,1))-diag(ones(N-1,1),-1);
B = zeros(N,1);

for i=1:N
    s = (tGLL(i)+tGLL(i+1))/2+(tGLL(i+1)-tGLL(i))/2*sG;
    h = LagrangeVal(s,N,1);
    
    B(i,1) = phi0*lambda*dt/2*difft(i)/2*sum(wG.*h(1,:));

    for k=1:N
        A(i,k) = A(i,k)-lambda*dt/2*difft(i)/2*sum(wG.*h(k+1,:));
    end
end
B(1,1) = B(1,1)+phi0;
    
phi = A\B;


pphi = [phi0 ; phi]'*hh;

plot(tt,pphi)


t = (tbegin+tend)/2+(tend-tbegin)/2*tGLL;

plot(t,[phi0 phi'],'o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error = sqrt( (phi_ex-pphi).^2);

subplot(1,2,2)
semilogy(tt,error,'-')
hold on
semilogy(tt(n),error(n),'or')
