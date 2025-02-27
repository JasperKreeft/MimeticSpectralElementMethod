

% 10 el T=10 dt=5e-2
% p-conv

N_case1 = [1 2 3 4 6 8]';
L2_case1 = [ 1.5807
             4.324e-1
             1.6821e-1
             6.0485e-2
             5.9097e-2
             5.8756e-2 ];

% h-conv N=2
% El_case2 = [ 10 20 30 40 50 ]';
El_case2 = [ 10 20 30 40 ]';
L2_case2 = [ 4.324e-1
             2.415e-1
             1.274e-1
             9.275e-2 ];
%              8.9819e-2 ];

% h-conv N=4
El_case3 = El_case2;
L2_case3 = [ 6.0485e-2
             7.7111e-2
             9.3003e-2
             1.0535e-1 ];
%              1.1680e-1 ];

% dt=5e-4
% p-conv
N_case4 = [2 3 4 6 8];
L2_case4 = [ 4.4329e-1
             1.3908e-1
             9.7753e-3
             6.2911e-3
             1.4948e-3 ];

% h-conv N=2
El_case5 = El_case2;
L2_case5 = [ 4.4329e-1
             2.7698e-1
             1.7982e-1
             1.2714e-1 ];

% h-conv N=4
El_case6 = El_case2;
L2_case6 = [ 9.7753e-3
             4.2155e-3
             1.9490e-3
             1.0367e-3 ];

% t-conv el=10 N=2
dt_case7 = [ 5e-2 5e-3 5e-4 ];
L2_case7 = [ 4.324e-1
             4.4329e-1
             4.4329e-1 ];
         
% t-conv el=10 N=4
dt_case8 = dt_case7;
L2_case8 = [ 6.0485e-2
             9.7753e-3
             9.7824e-3 ];

% 10 el T=10 dt=0.05
% p-conv
N_case9 = [2 3 4 6 8];
L2_case9 = [ 2.4150e-1
             8.3703e-2
             7.7111e-2
             7.7143e-2
             7.7110e-2 ];
         
% 10 el T=10 dt=5e-4
N_case10 = [ 2 4 8];
L2_case10 = [ 4.4329e-1
              9.7824e-3
              1.5100e-3 ];

% 10 el N=8 T=10
dt_case11 = dt_case7;
L2_case11 = [ 5.8756e-2
              1.4948e-3
              1.5100e-3 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
semilogy(N_case1,L2_case1,'-o')
hold on
semilogy(N_case4,L2_case4,'-sr')
semilogy(N_case10,L2_case10,'-^g')
legend('\Delta t=5e-2','\Delta t=5e-3','\Delta t=5e-4')
title('el=10 T=10')
xlabel('N')

figure
loglog(El_case2,L2_case2,'-o')
hold on
loglog(El_case3,L2_case3,'-sr')
legend('N=2','N=4')
title('T=10 \Delta t=5e-2')
xlabel('nr elements')

figure
semilogy(N_case1,L2_case1,'-o')
hold on
semilogy(N_case9,L2_case9,'-sr')
legend('10 el','20 el')
title('el=10 T=10 \Delta t=5e-2')
xlabel('N')

figure
loglog(El_case5,L2_case5,'-o')
hold on
loglog(El_case6,L2_case6,'-sr')
legend('N=2','N=4')
title('el=10 T=10 \Delta t=5e-4')
xlabel('nr elements')

figure
loglog(dt_case7,L2_case7,'-o')
hold on
loglog(dt_case8,L2_case8,'-sr')
legend('N=2','N=4')
title('T=10 el=10')
xlabel('\Delta t')