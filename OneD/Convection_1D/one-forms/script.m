clear all
close all
clc

N = 2;
H = 16%[ 2 4 8 16 32 ];
T = 100;  T = T-1e-10;
dt=5e-3;

for i=1:length(H)
    L2error_dt(i)=skew_one_1D_CN_ME_PerBC_funct(N,H(i),T,dt)
end

loglog(1./[16 32 64 128],[0.8034    0.3997    0.1742 0.0847])
hold on
loglog(1./[16 32 64 128],[0.0545    0.0223    0.0082 0.0029],'r')
loglog(1./[16 32 64 128],[0.0190    0.0059    0.0015 4.4998e-4],'g')
loglog(1./[8 16 32 64],[0.0139    0.0031    0.0011    0.0003],'c')
loglog(1./[4 8 16 32],[0.0335    0.0047    0.0010    0.0003],'m')

legend('N=1','N=2','N=3','N=4','N=5')
xlabel('\Delta x')
ylabel('L^2-error')
%%%%%%%%%%%%%%%%%%%%%%

% N = 20;
% H = 10;
% T = 1;  T = T-1e-10;
% dt=[5e-2 5e-3 5e-4];
% 
% for i=1:3
%     L2error_dt(i)=skew_one_1D_CN_ME_PerBC_funct(N,H,T,dt(i))
% end
% 
% N=40 H=10 9.0602e-003  6.0792e-005  8.5707e-006
% N=20 H=20 1.1982e-002  7.5492e-005  1.3687e-005

