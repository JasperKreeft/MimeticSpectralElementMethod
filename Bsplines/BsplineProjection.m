clear all
close all
clc

% Prange = 2:2:20;
Prange = 2;

% Hrange = [5 10 20 40 80 160 320 640 1280 2560];
Hrange = 6:4:120;
% Hrange = 9;

% Postproces grid and integration weights
[xgl,wgl] = GLLnodes(1000);

minX = 0;
maxX = 4;

x = (minX+maxX)/2+(maxX-minX)/2*xgl;

k=0;
error = zeros(1,max(length(Hrange),length(Prange)));
for P=Prange
for H=Hrange
k = k+1


%% 
X = linspace(0,4,H);
Jac = (maxX-minX)/2;

% Greville abscissa
Gr = GrevilleAbscissa(X,P);


% Function and its derivatives

% F(1,:) = cos(pi/2*Gr);
% F = zeros(P+1,length(Gr));
% for d=0:P
% F(d+1,:) = (-1)^ceil(d/2)*(pi/2)^d*(even(d)*cos(pi/2*Gr)+odd(d)*sin(pi/2*Gr));
% end

% % step function
F = zeros(P+1,length(Gr));
F(1,:) = (Gr>2);

% Reduction of function using dual functionals
L = reductionBspline(F,X,P);


%%

[B,E] = Bspline(x,X,P);

fh = L*B;

figure
% plot(x,cos(pi/2*x),'g')
plot(x,x>2.2,'g')
hold on
plot(Gr,L,'-sr')
plot(x,fh)


% fex = cos(pi/2*x);
fex = (x>2.2);

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