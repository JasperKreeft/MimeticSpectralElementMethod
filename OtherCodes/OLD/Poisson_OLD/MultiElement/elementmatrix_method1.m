function He = elementmatrix_method1(N,xiGLL,xiEG,wG,wGLL,XibGLLGLL,EtabGLLGLL)

global cc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gridtype = 'sinecurve'
dXdXibGLLGLL  = 1+pi*cc*cos(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
dXdEtabGLLGLL = pi*cc*sin(pi*XibGLLGLL).*cos(pi*EtabGLLGLL);
dYdXibGLLGLL  = pi*cc*cos(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
dYdEtabGLLGLL = 1+pi*cc*sin(pi*XibGLLGLL).*cos(pi*EtabGLLGLL);

dXibdXiGLLGLL   = (XibGLLGLL(N+1,1)-XibGLLGLL(1,1))/2*ones(N+1);
dXibdEtaGLLGLL  = zeros(N+1);
dEtabdXiGLLGLL  = zeros(N+1);
dEtabdEtaGLLGLL = (EtabGLLGLL(1,N+1)-EtabGLLGLL(1,1))/2*ones(N+1);

dXdXiGLLGLL  = dXdXibGLLGLL.*dXibdXiGLLGLL +dXdEtabGLLGLL.*dEtabdXiGLLGLL ;
dXdEtaGLLGLL = dXdXibGLLGLL.*dXibdEtaGLLGLL+dXdEtabGLLGLL.*dEtabdEtaGLLGLL;
dYdXiGLLGLL  = dYdXibGLLGLL.*dXibdXiGLLGLL +dYdEtabGLLGLL.*dEtabdXiGLLGLL ;
dYdEtaGLLGLL = dYdXibGLLGLL.*dXibdEtaGLLGLL+dYdEtabGLLGLL.*dEtabdEtaGLLGLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JGLLGLL =  dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = reshape(JGLLGLL,1,(N+1)^2);

Jacobian = spdiags(kron(j,ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

qinv11 = kron(reshape(( dXdXiGLLGLL./JGLLGLL),1,(N+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEtaGLLGLL./JGLLGLL),1,(N+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEtaGLLGLL./JGLLGLL),1,(N+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXiGLLGLL./JGLLGLL),1,(N+1)^2),[1 0])';

Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h  ,e  ] = MimeticpolyVal(xiGLL,N,1);
[hw ,ew ] = MimeticpolyVal(xiGLL,N,3);
[h_w,e_w] = MimeticpolyVal(xiEG ,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W1 = spdiags(kron(kron(wGLL,wGLL),ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

I1_GLLGLL = spalloc(2*(N+1)^2,2*N*(N+1),2*N*(N+1)^2);
I1_GLLGLL(1:2:2*(N+1)^2,1:N*(N+1)) = kron(e',speye(N+1));
for i=1:N
    I1_GLLGLL(2:2:2*(N+1)^2,(N+i-1)*(N+1)+(1:N+1)) = kron(speye(N+1),e(i,:)');
end

A = I1_GLLGLL'*Qinv'*(W1.*Jacobian)*Qinv*I1_GLLGLL;

I2_GLLG = kron(speye(2),kron(e_w(:,2:N+1),speye(N+1)));

W2 = spdiags([kron(wG,wGLL) kron(wG,wGLL)]',0,2*N*(N+1),2*N*(N+1));

I2_GGLL = kron(speye(2),kron(speye(N),ew));

B = I2_GLLG*W2*I2_GGLL';

He = inv(A)*B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%