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

L = M1*G/M0*G'*M1/(pi^2/4);

E = sort(eig(full(L),full(M1)));

F = E;
F(abs(F)<.1)=[];

plot(E,'.')

F(1:min(20,length(F)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
    rmpath('O:\MSEM\MSEM_codes\ThreeDim\Library');
else
    rmpath('/media/My Passport/MSEM/MSEM_codes/ThreeDim/Library');
end