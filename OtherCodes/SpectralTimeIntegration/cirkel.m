clear all
clf%ose all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)
x = linspace(-2,2,20);
x = x'*ones(size(x));
y=x';
quiver(x,y,-y,x,'c')
axis([-1 1 -1 1])
axis square
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 12;
tf = 2*pi;
x0 = 1;
y0 = 0;

t = GLLnodes(N);

[s0,w] = Gnodes(ceil((N+1)/2));


A = kron(eye(2),eye(N)-diag(ones(1,N-1),-1));
B = zeros(2*N,1);

for i=1:N

        s = (t(i)+t(i+1))/2+(t(i+1)-t(i))/2*s0;
        h = LagrangeVal(s,N,1);
        
        B(i,1)   = ((t(i+1)-t(i))/2*tf/2*sum(w.*h(1,:)));
        B(i+N,1) = B(i,1);
        
    for k=1:N
       
        A(i,k+N) = (t(i+1)-t(i))/2*tf/2*sum(w.*h(k+1,:));

        A(i+N,k) = -A(i,k+N);
    end
    
end

B = [-y0*ones(N,1) ; x0*ones(N,1)].*B;
B(1,1)   = B(1,1)   + x0;
B(N+1,1) = B(N+1,1) + y0;

gamma = A\B;

x = [x0 ; gamma(1:N)];
y = [y0 ; gamma(N+1:2*N)];

n = 1000;
tt = linspace(-1,1,n);
hh = LagrangeVal(tt,N,1);
xx = x'*hh;
yy = y'*hh;

subplot(1,2,1)
plot(xx,yy)
hold on
plot(x(1),y(1),'og')
plot(x(N+1),y(N+1),'or')
plot(x(2:N),y(2:N),'x')
axis([-2 2 -2 2])
axis square
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt_ex = linspace(0,tf,n);
exact = [cos(tt_ex) ; sin(tt_ex)];

error = sqrt( (exact(1,:)-xx).^2 + (exact(2,:)-yy).^2 );
% error = (error==0).*eps+error;

subplot(1,2,2)
semilogy(tt_ex,error,'-')
hold on
semilogy(tt_ex([1 end]),error([1 end]),'or')

