clear all
% close all
% clc

%% AFG2011 data Table 5.1

AFG_h = 1./[ 16 32 64 128 ];

AFG_U  = [ 3.26e-4 8.35e-5 2.10e-5 5.27e-6 ];
AFG_p  = [ 2.34e-3 8.05e-4 2.74e-4 9.39e-5 ];
AFG_w  = [ 2.70e-3 9.70e-4 3.47e-4 1.24e-4 ];
AFG_dw = [ 1.67e-1 1.24e-1 8.96e-2 6.42e-2 ];

%% KG 2012

KG_h = 1./[ 1 2 4 8 16 32 64 128 ];

KG_U  = [ 1.3023e-2 3.4750e-3 9.0356e-4 2.2699e-4 5.6808e-5 1.4206e-5 3.5516e-6 8.8793e-7 ];
KG_p  = [ 6.3905e-2 4.7041e-2 1.7859e-2 5.4311e-3 1.4932e-3 3.9121e-4 1.0010e-4 2.5318e-5 ];
KG_w  = [ 1.3041e-2 8.2541e-3 9.0512e-4 1.0280e-4 1.2445e-5 1.5424e-6 1.9238e-7 2.4035e-8 ];
KG_dw = [ 1.3452e-1 8.8122e-2 2.0897e-2 5.1304e-3 1.2773e-3 3.1903e-4 7.9741e-5 1.9934e-5 ];

KG_r_U  = [ 1.9060 1.9433 1.9930 1.9985 1.9996 1.9999 2.0000 ];
KG_r_p  = [ 0.4420 1.3973 1.7173 1.8628 1.9324 1.9664 1.9833 ];
KG_r_w  = [ 0.6599 3.1889 3.1383 3.0461 3.0124 3.0031 3.0008 ];
KG_r_dw = [ 0.6102 2.0762 2.0262 2.0060 2.0013 2.0003 2.0001 ];

KG_U_3  = [ 4.4535e-3 5.4116e-4 6.6826e-5 8.3274e-6 1.0401e-6 1.2999e-7 1.6248e-8  0 ];
KG_p_3  = [ 5.0912e-2 1.4807e-2 2.7680e-3 4.1493e-4 5.6544e-5 7.3720e-6 9.4086e-7  0 ];
KG_w_3  = [ 1.2919e-2 6.4247e-4 3.9421e-5 2.4609e-6 1.5379e-7 9.6121e-9 6.0076e-10 0 ];
KG_dw_3 = [ 1.2344e-1 1.2910e-2 1.5248e-3 1.8772e-4 2.3375e-5 2.9190e-6 3.6478e-7  0 ];

KG_r_U_3  = [ 3.0408 3.0176 3.0045 3.0011 3.0003 3.0001 3 ];
KG_r_p_3  = [ 1.7817 2.4194 2.7379 2.8754 2.9393 2.9700 3 ];
KG_r_w_3  = [ 4.3297 4.0266 4.0017 4.0001 4.0000 4.0000 4 ];
KG_r_dw_3 = [ 3.2573 3.0817 3.0220 3.0056 3.0014 3.0004 3 ];

%%

AX = [1e-3 1 1e-7 10];

figure
subplot(2,2,1)
loglog(AFG_h,AFG_U,'-o','markerface','b')
hold on
title('||u-u_h||')
loglog(KG_h(2:end),KG_U(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_U_3(2:end),'-og','markerface','g')
axis(AX)
legend('Raviart-Thomas, N=2','Mimetic Spectral, N=2','Mimetic Spectral, N=3',2)
set(legend,'box','off')
set(gca,'ytick',10.^(-12:2:2))
subplot(2,2,2)
loglog(AFG_h,AFG_p,'-o','markerface','b')
hold on
title('||p-p_h||')
loglog(KG_h(2:end),KG_p(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_p_3(2:end),'-og','markerface','g')
axis(AX)
set(gca,'ytick',10.^(-12:2:2))
subplot(2,2,3)
loglog(AFG_h,AFG_w,'-o','markerface','b')
hold on
title('||\omega-\omega_h||')
loglog(KG_h(2:end),KG_w(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_w_3(2:end),'-og','markerface','g')
axis(AX)
set(gca,'ytick',10.^(-12:2:2))
subplot(2,2,4)
loglog(AFG_h,AFG_dw,'-o','markerface','b')
hold on
title('||curl(\omega-\omega_h)||')
loglog(KG_h(2:end),KG_dw(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_dw_3(2:end),'-og','markerface','g')
axis(AX)
set(gca,'ytick',10.^(-12:2:2))