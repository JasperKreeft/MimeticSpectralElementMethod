clear all
close all
clc

%% AFG2011 data Table 2.1

AFG_h = 1./[ 16 32 64 128 ];

AFG_U  = [ 2.14e-3 5.37e-4 1.34e-4 3.36e-5 ];
AFG_dU = [ 1.17e-2 2.93e-3 7.33e-4 1.83e-4 ];
AFG_s  = [ 2.16e-4 2.70e-5 3.37e-6 4.16e-7 ];
AFG_ds = [ 2.63e-2 6.60e-3 1.65e-3 4.14e-4 ];

%% KG2012

KG_h = 1./[ 1 2 4 8 16 32 64 128 ];

KG_U  = [ 5.9664e-1 1.0665e-1 2.5957e-2 6.4417e-3 1.6074e-3 4.0166e-4 1.0040e-4 2.5100e-5 ];
KG_dU = [ 2.8419e+0 6.1409e-1 1.5354e-1 3.8326e-2 9.5769e-3 2.3939e-3 5.9846e-4 1.4961e-4 ];
KG_s  = [ 5.6723e-1 6.2254e-2 6.5422e-3 7.8395e-4 9.7018e-5 1.2098e-5 1.5113e-6 1.8888e-7 ];
KG_ds = [ 3.3202e+0 7.7839e-1 1.7100e-1 4.0808e-2 1.0072e-2 2.5096e-3 6.2687e-4 1.5669e-4 ];

KG_r_U  = [ 2.4840 2.0386 2.0106 2.0027 2.0007 2.0002 2.0000 ];
KG_r_dU = [ 2.2103 1.9998 2.0022 2.0007 2.0002 2.0000 2.0000 ];
KG_r_s  = [ 3.1877 3.2503 3.0609 3.0144 3.0035 3.0009 3.0002 ];
KG_r_ds = [ 2.0927 2.1865 2.0671 2.0186 2.0048 2.0012 2.0003 ];

KG_U_3  = [ 5.7011e-2 1.3636e-2 1.7060e-3 2.1326e-4 2.6657e-5 3.3321e-6 4.1652e-7 5.2065e-8 ];
KG_dU_3 = [ 3.1275e-1 8.0122e-2 1.0134e-2 1.2701e-3 1.5886e-4 1.9861e-5 2.4827e-6 3.1034e-7 ];
KG_s_3  = [ 4.0266e-2 4.9624e-3 2.8890e-4 1.7672e-5 1.0983e-6 6.8547e-8 4.2827e-9 2.6766e-10 ];
KG_ds_3 = [ 5.1997e-1 8.9119e-2 1.0789e-2 1.3357e-3 1.6654e-4 2.0804e-5 2.6001e-6 3.2500e-7 ];

KG_r_U_3  = [ 2.0638 2.9987 2.9999 3.0000 3.0000 3.0000 3.0000 ];
KG_r_dU_3 = [ 1.9647 2.9831 2.9962 2.9991 2.9998 2.9999 3.0000 ];
KG_r_s_3  = [ 3.0205 4.1024 4.0311 4.0081 4.0020 4.0005 4.0000 ];
KG_r_ds_3 = [ 2.5446 3.0462 3.0139 3.0036 3.0009 3.0002 3.0001 ];


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
subplot(2,2,2)
loglog(AFG_h,AFG_dU,'-o','markerface','b')
hold on
title('||div(u-u_h)||')
loglog(KG_h(2:end),KG_dU(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_dU_3(2:end),'-og','markerface','g')
axis(AX)
subplot(2,2,3)
loglog(AFG_h,AFG_s,'-o','markerface','b')
hold on
title('||\sigma-\sigma_h||')
loglog(KG_h(2:end),KG_s(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_s_3(2:end),'-og','markerface','g')
axis(AX)
subplot(2,2,4)
loglog(AFG_h,AFG_ds,'-o','markerface','b')
hold on
title('||curl(\sigma-\sigma_h)||')
loglog(KG_h(2:end),KG_ds(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_ds_3(2:end),'-og','markerface','g')
axis(AX)

%% AFG2011 data Table 2.2

AFG_h = 1./[ 16 32 64 128 ];

AFG_U  = [ 1.22e-3 3.05e-4 7.63e-5 1.91e-5 ];
AFG_dU = [ 1.55e-2 5.33e-3 1.85e-3 6.49e-4 ];
AFG_s  = [ 1.90e-2 6.36e-3 2.18e-3 7.58e-4 ];
AFG_ds = [ 2.53e-0 1.68e-0 1.14e-0 7.89e-1 ];




%% KG2012

KG_h = 1./[ 1 2 4 8 16 32 64 128 ];

KG_U  = [ 3.0897e-1 6.7449e-2 1.6417e-2 4.0741e-3 1.0166e-3 2.5403e-4 6.3500e-5 1.5875e-5 ];
KG_dU = [ 9.9855e-1 2.8949e-1 7.2380e-2 1.8067e-2 4.5146e-3 1.1285e-3 2.8212e-4 7.0528e-5 ];
KG_s  = [ 2.8349e-1 7.2041e-2 8.8157e-3 1.0960e-3 1.3682e-4 1.7096e-5 2.1369e-6 2.6710e-7 ];
KG_ds = [ 3.6070e+0 9.2554e-1 2.2834e-1 5.6818e-2 1.4187e-2 3.5455e-3 8.8631e-4 2.2157e-4 ];

KG_r_U  = [ 2.1956 2.0386 2.0106 2.0027 2.0007 2.0002 2.0000 ];
KG_r_dU = [ 1.7863 1.9998 2.0022 2.0007 2.0002 2.0000 2.0000 ];
KG_r_s  = [ 1.9764 3.0307 3.0078 3.0019 3.0005 3.0001 3.0000 ];
KG_r_ds = [ 1.9624 2.0191 2.0067 2.0018 2.0005 2.0001 2.0000 ];

KG_U_3  = [ 3.9989e-2 8.6240e-3 1.0789e-3 1.3488e-4 1.6859e-5 2.1074e-6 2.6343e-7 3.2929e-8 ];
KG_dU_3 = [ 2.8294e-1 3.7770e-2 4.7770e-3 5.9872e-4 7.4889e-5 9.3626e-6 1.1704e-6 1.4630e-7 ];
KG_s_3  = [ 1.0110e-1 6.2854e-3 3.9586e-4 2.4789e-5 1.5501e-6 9.6890e-8 6.0558e-9 3.7850e-10 ];
KG_ds_3 = [ 9.1874e-1 1.1958e-1 1.5035e-2 1.8818e-3 2.3530e-4 2.9414e-5 3.6769e-6 4.5961e-7 ];

KG_r_U_3  = [ 2.2132 2.9987 2.9999 3.0000 3.0000 3.0000 3.0000 ];
KG_r_dU_3 = [ 2.9052 2.9831 2.9962 2.9991 2.9998 2.9999 3.0000 ];
KG_r_s_3  = [ 4.0077 3.9889 3.9972 3.9993 3.9998 4.0000 4.0000 ];
KG_r_ds_3 = [ 2.9416 2.9916 2.9981 2.9996 2.9999 3.0000 3.0000 ];

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
subplot(2,2,2)
loglog(AFG_h,AFG_dU,'-o','markerface','b')
hold on
title('||div(u-u_h)||')
loglog(KG_h(2:end),KG_dU(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_dU_3(2:end),'-og','markerface','g')
axis(AX)
subplot(2,2,3)
loglog(AFG_h,AFG_s,'-o','markerface','b')
hold on
title('||\sigma-\sigma_h||')
loglog(KG_h(2:end),KG_s(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_s_3(2:end),'-og','markerface','g')
axis(AX)
subplot(2,2,4)
loglog(AFG_h,AFG_ds,'-o','markerface','b')
hold on
title('||curl(\sigma-\sigma_h)||')
loglog(KG_h(2:end),KG_ds(2:end),'-or','markerface','r')
loglog(KG_h(2:end),KG_ds_3(2:end),'-og','markerface','g')
axis(AX)



