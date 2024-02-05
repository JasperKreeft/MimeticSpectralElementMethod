function [A,B,C] = som_d1p0(N,xiGLL,xiEG,wG,wGLL,NGp)
% clear all; close all; clc
% N=2;
% [xiG,wG]=Gnodes(N);
% [xiGLL,wGLL]=GLLnodes(N);
% xiEG = [-1 xiG 1];

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);

e    = EdgeVal(dhdxi );
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wgl = spdiags(kron(wGLL,wGLL)',0,(N+1)^2,(N+1)^2);

A = Wgl;

Wgl_bc = ??????
Iggl_bc = ?????????
B=?????????????




Iglg = [ kron(speye(N+1),e_w(:,2:N+1))                zeros(N*(N+1))
                        zeros(N*(N+1)) kron(e_w(:,2:N+1),speye(N+1)) ];

Iggl = [ kron(speye(N),ew)    zeros(N*(N+1))
         zeros(N*(N+1))    kron(ew,speye(N)) ];

Wglg = spdiags([kron(wG,wGLL) kron(wGLL,wG)]',0,2*N*(N+1),2*N*(N+1));

P = [ zeros(N*(N+1)) -speye(N*(N+1))
      speye(N*(N+1))  zeros(N*(N+1)) ];

C = NGp*Iglg*P*Wglg*Iggl';