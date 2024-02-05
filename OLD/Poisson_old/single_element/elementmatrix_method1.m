function [A,B,Dp,Gd,F] = elementmatrix_method1(N,xiGLL,xiEG,wG,wGLL,XiGLLGLL,EtaGLLGLL)

global gridtype c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch gridtype
    case 'standard'

        dXdXiGLLGLL  = ones(N+1);  dXdEtaGLLGLL = zeros(N+1);
        dYdXiGLLGLL  = zeros(N+1); dYdEtaGLLGLL = ones(N+1);
        
    case 'sinecurve'
       
        dXdXiGLLGLL  = 1+pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
        dXdEtaGLLGLL = pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
        dYdXiGLLGLL  = pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
        dYdEtaGLLGLL = 1+pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JGLLGLL = dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = reshape(JGLLGLL',1,(N+1)^2);

Jacobian = sparse(diag([j j]));

Qinv = sparse([ diag(reshape(( dXdXiGLLGLL./JGLLGLL),1,(N+1)^2)) diag(reshape(( dXdEtaGLLGLL./JGLLGLL),1,(N+1)^2))
                diag(reshape(( dYdXiGLLGLL./JGLLGLL),1,(N+1)^2)) diag(reshape(( dYdEtaGLLGLL./JGLLGLL),1,(N+1)^2)) ]); %zeros((N+1)^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);

e    = EdgeVal(dhdxi );
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);

W1 = spdiags([kron(wGLL,wGLL) kron(wGLL,wGLL)]',0,2*(N+1)^2,2*(N+1)^2); % is [kron(wGLL_q,wGLL_p) kron(wGLL_p,wGLL_q)]

I1_GLLGLL_xi = kron(speye(N+1),e');
I1_GLLGLL_eta = kron(e',speye(N+1));

I1_GLLGLL = [ I1_GLLGLL_eta zeros((N+1)^2,N*(N+1)) ; zeros((N+1)^2,N*(N+1)) I1_GLLGLL_xi ];

A = I1_GLLGLL'*Qinv'*(W1.*Jacobian)*Qinv*I1_GLLGLL;

I2_GLLG_xi  = kron(speye(N+1),e_w(:,2:N+1));
I2_GLLG_eta = kron(e_w(:,2:N+1),speye(N+1));

I2_GLLG = [ I2_GLLG_eta zeros(N*(N+1)) ; zeros(N*(N+1)) I2_GLLG_xi ];

W2 = spdiags([kron(wG,wGLL) kron(wGLL,wG)]',0,2*N*(N+1),2*N*(N+1));

I2_GGLL_xi  = kron(speye(N),ew');
I2_GGLL_eta = kron(ew',speye(N));

I2_GGLL = [ I2_GGLL_xi zeros(N*(N+1)) ; zeros(N*(N+1)) I2_GGLL_eta ];

B = I2_GLLG*W2*I2_GGLL;

[Dp,Gd] = topology(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = reshape(force(xiGLL,xiGLL),N*N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%