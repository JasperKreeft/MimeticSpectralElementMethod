clear all
close all
clc

if ispc
    path(path,'O:\MSEM\MSEM_codes\ThreeDim\Library');
else
    path(path,'/media/My Passport/MSEM/MSEM_codes/ThreeDim/Library');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N w e

N = 5;

[xi,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(xi,N,1);

G = grad_out_3D(N);

M0 = innerproduct_3D(0);
M1 = innerproduct_3D(1);

L = G'*M1*G/(pi^2/4);

E = sort(eig(full(L),full(M0)));

E(1:min(20,length(E)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    rmpath('O:\MSEM\MSEM_codes\ThreeDim\Library');
else
    rmpath('/media/My Passport/MSEM/MSEM_codes/ThreeDim/Library');
end