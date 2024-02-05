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

N = 6;

[xi,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(xi,N,1);

C = curl_out_3D(N);
D = div_out_3D(N);

M1 = innerproduct_3D(1);
M2 = innerproduct_3D(2);
M3 = innerproduct_3D(3);

L = M2*C/M1*C'*M2+D'*M3*D;

E = sort(eig(full(L),full(M2*(pi^2/4))));

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