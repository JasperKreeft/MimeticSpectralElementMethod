function [A,B,R] = elementmatrix_singlegrid(Xib,Etab)

global N xiGLL wGLL cc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% gridtype = 'sinecurve'
dXdXib  = 1+pi*cc*cos(pi*Xib).*sin(pi*Etab);
dXdEtab = pi*cc*sin(pi*Xib).*cos(pi*Etab);
dYdXib  = pi*cc*cos(pi*Xib).*sin(pi*Etab);
dYdEtab = 1+pi*cc*sin(pi*Xib).*cos(pi*Etab);

dXibdXi   = (Xib(N+1,1)-Xib(1,1))/2*ones(N+1);
dXibdEta  = zeros(N+1);
dEtabdXi  = zeros(N+1);
dEtabdEta = (Etab(1,N+1)-Etab(1,1))/2*ones(N+1);

dXdXi  = dXdXib.*dXibdXi +dXdEtab.*dEtabdXi ;
dXdEta = dXdXib.*dXibdEta+dXdEtab.*dEtabdEta;
dYdXi  = dYdXib.*dXibdXi +dYdEtab.*dEtabdXi ;
dYdEta = dYdXib.*dXibdEta+dYdEtab.*dEtabdEta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = dXdXi.*dYdEta-dXdEta.*dYdXi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = reshape(J,1,(N+1)^2);
% keyboard
Jacobian = spdiags([j j]',0,2*(N+1)^2,2*(N+1)^2);

qinv11 = kron(reshape(( dXdXi./J),1,(N+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEta./J),1,(N+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEta./J),1,(N+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXi./J),1,(N+1)^2),[1 0])';

Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xiGLL,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W1 = spdiags(kron(kron(wGLL,wGLL),ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

I1 = spalloc(2*(N+1)^2,2*N*(N+1),2*N*(N+1)^2);
I1(1:2:2*(N+1)^2,1:N*(N+1)) = kron(e',speye(N+1));
for i=1:N
    I1(2:2:2*(N+1)^2,(N+i-1)*(N+1)+(1:N+1)) = kron(speye(N+1),e(i,:)');
end

% inner-product matrix for 1-forms ( or (n-1)-forms with n=2 )
A = I1'*Qinv'*(W1.*Jacobian)*Qinv*I1;

% inner-product matrix for 2-forms ( or n-forms with n=2 )
B = zeros(N^2);
for k=1:N
    for l=1:N
        kl = k+(l-1)*N;
        for i=1:N
            for j=1:N
                ij = i+(j-1)*N;
                B(kl,ij) = B(kl,ij) + (wGLL.*e(i,:).*e(k,:))*J*(wGLL.*e(j,:).*e(l,:))';
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions

R = zeros(2*N*(N+1),N^2);
for i=1:N
    for j=1:N
        ij = i+(j-1)*N;
        for l=1:N
            int = sum(wGLL.*e(j,:).*e(l,:));
            kl = 1+(l-1)*(N+1);
            R(kl,ij) = R(kl,ij)-e(i,1)*int;
            kl = l*(N+1);
            R(kl,ij) =R(kl,ij) +e(i,N+1)*int;
        end
        for k=1:N
            int = sum(wGLL.*e(i,:).*e(k,:));
            kl = N*(N+1)+1+(k-1)*(N+1);
            R(kl,ij) = R(kl,ij)-e(j,1)*int;
            kl = N*(N+1)+k*(N+1);
            R(kl,ij) = R(kl,ij)+e(j,N+1)*int;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%