clear all
close all
clc

N=1;

[xGLL,wGLL] = GLLnodes(N);
[xG,wG] = Gnodes(N);
xEG = [-1 xG 1];
wEG = [0 wG 0];
    
[h   ,dhdxi  ] = LagrangeVal(xGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xEG,N,1);
[hw_w,dhwdxiw] = LagrangeVal(xEG,N,3);

e    = EdgeVal(dhdxi);
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);
ew_w = EdgeVal(dhwdxiw);

in_qq = diag(wGLL);
ex_uq = (wGLL'*ones(1,N+1)).*ew';
ex_pq = -dhdxiw.*(ones(N+1,1)*wEG);
ex_pq(1,1) = ex_pq(1,1)-1;
ex_pq(N+1,N+2) = ex_pq(N+1,N+2)+1;

E10 = [diag(-ones(1,N+1)) zeros(N+1,1)]+[zeros(N+1,1) diag(ones(1,N+1))];

ex_uq_p = ex_uq*E10;

figure
subplot(2,2,1)
surf(in_qq)
% view([0 0 1])
subplot(2,2,2)
surf(ex_uq_p)
% view([0 0 1])
subplot(2,2,3)
surf(ex_pq)
subplot(2,2,4)
surf(abs(ex_uq_p-ex_pq))