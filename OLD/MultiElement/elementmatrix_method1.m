function [Dp,H,Gd] = elementmatrix_method1(N,xiGLL,xiEG,wG,wGLL,numRows,numColumns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);

e    = EdgeVal(dhdxi );
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);

J = 1/(numRows*numColumns);

dxdxi  = 1/numColumns;
dydeta = 1/numRows;

Jacobian = sparse(diag(J*ones(1,2*(N+1)^2)));

Qinv =??????;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W1 = sparse(diag([kron(wGLL,wGLL) kron(wGLL,wGLL)])); % is [kron(wGLL_q,wGLL_p) kron(wGLL_p,wGLL_q)]

I1_GLLGLL_eta = kron(e',eye(N+1));
I1_GLLGLL_xi = kron(eye(N+1),e');

I1_GLLGLL = [ I1_GLLGLL_eta zeros((N+1)^2,N*(N+1)) ; zeros((N+1)^2,N*(N+1)) I1_GLLGLL_xi ];

A = I1_GLLGLL'*Qinv'*(W1.*Jacobian)*Qinv*I1_GLLGLL;

I2_GLLG_eta = kron(e_w(:,2:N+1),eye(N+1));
I2_GLLG_xi = kron(eye(N+1),e_w(:,2:N+1));

I2_GLLG = [ I2_GLLG_eta zeros(N*(N+1)) ; zeros(N*(N+1)) I2_GLLG_xi ];

W2 = sparse(diag([kron(wG,wGLL) kron(wGLL,wG)]));

I2_GGLL_xi  = kron(eye(N),ew');
I2_GGLL_eta = kron(ew',eye(N));

I2_GGLL = [ I2_GGLL_xi zeros(N*(N+1)) ; zeros(N*(N+1)) I2_GGLL_eta ];

B = I2_GLLG*W2*I2_GGLL;

H = inv(A)*B;

[Dp,Gd] = topology(N);

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%