clear all
close all
clc

Prange = 2:2:20;
% Prange = 2;

% Hrange = [5 10 20 40 80 160 320 640 1280 2560];
% Hrange = 6:4:120;
Hrange = 5;

k=0;
error = zeros(1,max(length(Hrange),length(Prange)));
for P=Prange
for H=Hrange
k = k+1
    
X = linspace(0,4,H);
Xi = [ X(1)*ones(1,P) X X(end)*ones(1,P) ];

%% Greville abscissa
Gr = zeros(1,length(X)+P-1);
for j=1:length(X)+P-1
    Gr(j) = 1/P*sum(Xi(j+1:j+P));
end

% F(1,:) = cos(pi/2*Gr);
F = zeros(P+1,length(Gr));
for d=0:P
F(d+1,:) = (-1)^ceil(d/2)*(pi/2)^d*(even(d)*cos(pi/2*Gr)+odd(d)*sin(pi/2*Gr));
end

% % step function
% F = zeros(P+1,length(Gr));
% F(1,:) = (Gr>2);



%% 
L = zeros(1,length(X)+P-1);
for j=1:length(X)+P-1

xi_j_jp = Xi(j+1:j+P);
n = length(xi_j_jp)+1;
phi = zeros(P+1,n);
phi(1,:) = poly(xi_j_jp);
for a = 1:P
    phi(a+1,:) = [0 (n-1:-1:1).*phi(a,1:n-1)];
end

phiGr = zeros(P+1,1);
for a=1:P+1
    phiGr(a) = phi(a,:)*(Gr(j).^(size(phi,2)-1:-1:0))';
end


for d=0:P
    L(j) = L(j) + (-1)^d/factorial(P)*phiGr(P-d+1)*F(d+1,j);
end

end

%%

% x = linspace(min(X),max(X),1000);
[xgl,wgl] = GLLnodes(100);
x = (min(X)+max(X))/2+(max(X)-min(X))/2*xgl;

[B,E] = Bspline(x,X,P);

fh = L*B;

% figure
% plot(x,cos(pi/2*x),'g')
% hold on
% plot(Gr,L,'-sr')
% plot(x,fh)


fex = cos(pi/2*x);
% fex = (x>2);

Jac = 4/2;
error(k) = sqrt(sum((fex-fh).^2.*wgl*Jac));

end
end
%%
if k>1
figure
xaxis = 4./(Hrange-1);
if length(Hrange)>1
loglog(xaxis,error,'-ok','markerface','k')
ylim([1e-10 1e-1])
xlabel('h')
set(gca,'ytick',10.^(-10:2:0))
rate = (log(error(k))-log(error(3)))/(log(xaxis(end))-log(xaxis(3)));
text(5e-2,1e-7,['rate = ' num2str(roundn(rate,-1))])
end
if length(Prange)>1
semilogy(Prange,error,'-ok','markerface','k')
ylim([1e-10 1e-1])
set(gca,'ytick',10.^(-10:2:0))
xlabel('p')
end
ylabel('L^2 error')
end