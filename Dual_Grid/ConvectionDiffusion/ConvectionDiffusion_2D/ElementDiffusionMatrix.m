function [Diff_1,Diff_0] = ElementSupportOperatorMethod()

global N N2
global nodes_in_element edges_in_element
global wG wGLL
global e e_w

int_point = (N+1)^2; % integration points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W1 = spdiags(kron(kron(wGLL,wGLL),ones(1,2))',0,2*int_point,2*int_point);

I1_GLLGLL = spalloc(2*int_point,edges_in_element,2*N*(N+1)^2);
I1_GLLGLL(1:2:2*int_point,1:edges_in_element) = kron(e',speye(N+1));
for i=1:N
    I1_GLLGLL(2:2:2*int_point,(N+i-1)*(N+1)+(1:N+1)) = kron(speye(N+1),e(i,:)');
end

Diff_1 = I1_GLLGLL'*Qinv'*(W1.*Jacobian)*Qinv*I1_GLLGLL;

Ii_GLLG_eta = kron(speye(N),e_w(:,2:N+1));
Ii_GLLG_xi  = kron(e_w(:,2:N+1),speye(N));

% Boundary condition weights !!!!!!!!!!!!!!
Ib = kron(speye(2),kron(e_w(:,2:N+1),speye(2)));

I2 = spalloc(nodes_in_element,nodes_in_element,N2*nodes_in_element);
I2(1:N2,1:N2)             = Ii_GLLG_xi*Ii_GLLG_eta;
I2(N2+(1:4*N),N2+(1:4*N)) = Ib;

Wi = kron(wG,wG);
Wb = kron([1 1],kron(wG,[1 1]));

W2 = spdiags([Wi Wb]',0,nodes_in_element,nodes_in_element);

Diff_0 = I2*W2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%