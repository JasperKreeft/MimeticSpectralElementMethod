function [A,B] = hodge_p1d1(N,xiGLL,xiEG,wG,wGLL)
% clear all; close all; clc
% N=2;
% [xiG,wG]=Gnodes(N);
% [xiGLL,wGLL]=GLLnodes(N);
% xiEG = [-1 xiG 1];

[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);
[hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);

ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);
ew_w = EdgeVal(dhwdxiw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wgg = spdiags([kron(wG,wG) kron(wG,wG)]',0,2*N^2,2*N^2);

Igg = [     kron(speye(N),ew_w(:,2:N+1)) zeros(N*(N+1),N*N)
          zeros(N*(N+1),N*N)     kron(ew_w(:,2:N+1),speye(N)) ];

A = Igg*Wgg*Igg';

Iglg = [ kron(speye(N+1),e_w(:,2:N+1))                zeros(N*(N+1))
                        zeros(N*(N+1)) kron(e_w(:,2:N+1),speye(N+1)) ];

Iggl = [ kron(speye(N),ew)    zeros(N*(N+1))
         zeros(N*(N+1))    kron(ew,speye(N)) ];

Wglg = spdiags([kron(wG,wGLL) kron(wGLL,wG)]',0,2*N*(N+1),2*N*(N+1));

P = [ zeros(N*(N+1)) -speye(N*(N+1))
      speye(N*(N+1))  zeros(N*(N+1)) ];

B = Iggl*Wglg*P*Iglg';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%