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

N = 8;

[xi,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(xi,N,1);

D = div_out_3D(N);

M2 = innerproduct_3D(2);
M3 = innerproduct_3D(3);

L = M3*D/M2*D'*M3/(pi^2/4);

E = sort(eig(full(L),full(M3)));

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