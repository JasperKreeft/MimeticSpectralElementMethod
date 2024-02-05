%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Old code. Predicesor of fv_sem_v2.m                                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

mmax = 1;
color = 'brgcmky';
leg = zeros(mmax,4);
for m=1:mmax

ZZ = 3:10;%10:20;
errorL1 = zeros(size(ZZ));
errorL2 = zeros(size(ZZ));
for Z=ZZ

N = Z

gridchoice = 2; % 1 uniform, 2 non-uniform

if gridchoice == 1
    xp = linspace(-1,1,N+1);
    xd = linspace((xp(1)+xp(2))/2,(xp(N)+xp(N+1))/2,N);
    xd_ex = [-1 xd 1];
    [LaPoly,dLaPoly,Pweights] = LagrangePoly(xp);        % primal polynomials
    [LaPolyw,dLaPolyw,Dexweights] = LagrangePoly(xd_ex); % dual polynomials (w = wiggle)
    Dweights = Dexweights(2:N+1);
elseif gridchoice == 2
    [xp,Pweights] = gausslobattolegendre(N);
    [xd,Dweights] = gauss(N);
    [xd_ex,Dexweights] = extendedgauss(N);
    [LaPoly,dLaPoly] = LagrangePoly(xp);      % primal polynomials
    [LaPolyw,dLaPolyw] = LagrangePoly(xd_ex); % dual polynomials (w = wiggle)
end

yp = xp; yd = xd; yd_ex = xd_ex;
dxd_ex = diff(xd_ex); dyd_ex = dxd_ex;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dh_i/dx (x_j)
dhdx = zeros(N+1);
for i=1:N+1
    dhdx(i,:) = polyval_J(dLaPoly(i,:),xp);
end

% dhw_i/dx (x_j)
dhwdx = zeros(N+2,N+1);
for i=1:N+2
    dhwdx(i,:) = polyval_J(dLaPolyw(i,:),xp);
end

% dh_i/dx (xw_j)
dhdxw = zeros(N+1,N+2);
for i=1:N+1
    dhdxw(i,:) = polyval_J(dLaPoly(i,:),xd_ex);
end

% dhw_i/dx (xw_j)
dhwdxw = zeros(N+2);
for i=1:N+2
    dhwdxw(i,:) = polyval_J(dLaPolyw(i,:),xd_ex);
end


% e_i(x_j) = - sum_{k=0}^{i-1} dh_k/dx(x_j)
e = zeros(N,N+1);
for i=1:N
    for j=1:N+1
        e(i,j) = -sum(dhdx(1:i,j));
    end
end

% ew_i(x_j) = - sum_{k=0}^{i-1} dhw_k/dx(x_j)
ew = zeros(N+1);
for i=1:N+1
    for j=1:N+1
        ew(i,j) = -sum(dhwdx(1:i,j));
    end
end

% e_i(xw_j) = - sum_{k=0}^{i-1} dh_k/dx(xw_j)
e_w = zeros(N,N+2);
for i=1:N
    for j=1:N+2
        e_w(i,j) = -sum(dhdxw(1:i,j));
    end
end

% ew_i(xw_j) = - sum_{k=0}^{i-1} dhw_k/dx(xw_j)
ew_w = zeros(N+1,N+2);
for i=1:N+1
    for j=1:N+2
        ew_w(i,j) = -sum(dhwdxw(1:i,j));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Au = zeros(N*(N+1)); Bu = zeros(N*(N+1));
for l=1:N
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        % Au
        for j=1:N+2
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                Au(kl,ij) = ew(i,k)*e_w(l,j)*Pweights(k)*Dexweights(j);
            end
        end
        
        %Bu
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                Bu(kl,ij) = (i==k)*Pweights(i)*sum(Pweights.*e(j,:).*e(l,:));
            end
        end
        
    end
end


Av = zeros(N*(N+1)); Bv = zeros(N*(N+1));
for l=1:N+1
    for k=1:N
        kl = k+(l-1)*N;
        
        % Av
        for j=1:N+1
            for i=1:N+2
                ij  = i+(j-1)*(N+2);
                Av(kl,ij) = e_w(k,i)*ew(j,l)*Dexweights(i)*Pweights(l);
            end
        end

        % Bv
        for j=1:N+1
            for i=1:N
                ij  = i+(j-1)*N;
                Bv(kl,ij) = (j==l)*Pweights(j)*sum(Pweights.*e(i,:).*e(k,:));
            end
        end

    end
end


Au(:,[1:N+1 (N+1)*(N+1)+1:(N+1)*(N+2)])=[];
Av(:,[1:(N+2):(N+2)*(N+1) (N+2):(N+2):(N+2)*(N+1)]) = [];


H = [inv(Bu)*Au zeros(N*(N+1)); zeros(N*(N+1)) inv(Bv)*Av];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = zeros(N*N,1);
for j=1:N
    for i=1:N
        k = i+(j-1)*N;
        F(k) = 2*(cos(m*pi*xp(i+1))-cos(m*pi*xp(i)))*(cos(m*pi*yp(j+1))-cos(m*pi*yp(j)));
    end
end
phi_exact = zeros((N+2)*(N+2),1);
for j=1:N+2
    for i=1:N+2
        k = i+(j-1)*(N+2);
        phi_exact(k) = sin(m*pi*xd_ex(i))*sin(m*pi*yd_ex(j));
    end
end

PHI_EXACT = zeros(N+2);
for i=1:N+2
    for j=1:N+2
        k = i+(j-1)*(N+2);
        PHI_EXACT(i,j) = phi_exact(k);
    end
end

Dp = topology(N);
Gd = Dp';

A = Dp*H*Gd;


% Additing boundary conditions
% for i=1:N
%     for j=1:N
%         k = i+(j-1)*N;
%         if i==1
%             F(k,1) = F(k,1) + 2*sin(m*pi*xd_ex(1))*sin(m*pi*yd_ex(j+1));
%         elseif i==N
%             F(k,1) = F(k,1) + 2*sin(m*pi*xd_ex(N+2))*sin(m*pi*yd_ex(j+1));
%         end
%         if j==1
%             F(k,1) = F(k,1) + 2*sin(m*pi*xd_ex(i+1))*sin(m*pi*yd_ex(1));
%         elseif j==N
%             F(k,1) = F(k,1) + 2*sin(m*pi*xd_ex(i+1))*sin(m*pi*yd_ex(N+2));
%         end
%     end
% end

phi_in = A\F;

PHI = zeros(N+2); l=0;
for i=1:N+2
    for j=1:N+2
        k = i+(j-1)*(N+2);
        if i==1 || i==N+2
            PHI(i,j) = phi_exact(k);
        elseif j==1 || j==N+2
            PHI(i,j) = phi_exact(k);
        else            
        l=l+1;
        PHI(i,j) = phi_in(l);
        end
    end
end

% surf(xd_ex,yd_ex,PHI)

%% Exact

% xx = linspace(-1,1,100);
% yy = linspace(0-1,1,100);
% 
% phi_ex = zeros(100);
% for i=1:100
%     for j=1:100
%         phi_ex(i,j) = sin(m*pi*xx(i))*sin(m*pi*yy(j));
%     end
% end
% 
% hold on
% 
% surf(xx,yy,phi_ex')

%% error

errorL1(N+1-ZZ(1)) = 1/N^2*sum(sum(abs(PHI-PHI_EXACT)));
errorL2(N+1-ZZ(1)) = sqrt(1/N^2*sum(sum((PHI-PHI_EXACT).^2)));
end

figure(1)
semilogy(ZZ,errorL1,color(m))
hold on
xlim([ZZ(1) N])
xlabel('N')
ylabel('L^1-error')

figure(2)
semilogy(ZZ,errorL2,color(m))
hold on
xlim([ZZ(1) N])
xlabel('N')
ylabel('L^2-error')

leg(m,:) = 'm = ';
end
figure(1)
legend([leg num2str((1:mmax)')])
figure(2)
legend([leg num2str((1:mmax)')])