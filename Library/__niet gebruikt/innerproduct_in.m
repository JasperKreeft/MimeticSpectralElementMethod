function M = innerproduct_in(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Inner-product matrices for k-forms, k=0,1,2                             %
% for a single spectral element                                           %
%                                                                         %
% M0 = innerproduct_in(0,J)                                               %
% M1 = innerproduct_in(1,J,Qinv)                                          %
% M2 = innerproduct_in(2,J)                                               %
%                                                                         %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 28-10-2010                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N w e

form = varargin{1};

if nargin<2
    warning('JK: Not enough input arguments')
    J = ones(N+1);
else
    J = varargin{2};
end



if size(J) == [ N+1 N+1 ]                                      %#ok<*BDSCA>
    J = reshape(J,(N+1)^2,1);
elseif size(J) == [ 1 (N+1)^2 ]
    J = J';
end


switch form

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Inner-product matrix for 0-forms (or (n-2)-forms with n=2)              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 0

W0 = kron(w,w)';
M = spdiags(J.*W0,0,(N+1)^2,(N+1)^2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Inner-product matrix for 1-forms (or (n-1)-forms with n=2)              %
% for a single spectral element                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 1

if form==1 && nargin==3
    Qinv = varargin{3};
else
    error('JK: Number of input arguments are not correct')
end

if size(Qinv) == [2*(N+1)^2 3]
    Qinv = spdiags(Qinv,-1:1,2*(N+1)^2,2*(N+1)^2);
end

W1  = kron(w,w)';
% W1 = kron(w,ones(1,N+1))';
JW1 = spdiags(kron(ones(2,1),J.*W1),0,2*(N+1)^2,2*(N+1)^2);

% I1 = spalloc(2*(N+1)^2,2*N*(N+1),2*N*(N+1)^2);
% I1(1:2:2*(N+1)^2,1:N*(N+1)) = kron(e',speye(N+1));
% for i=1:N
%     I1(2:2:2*(N+1)^2,(N+i-1)*(N+1)+(1:N+1)) = kron(speye(N+1),e(i,:)');
% end

I1 = zeros(2*(N+1)^2,2*N*(N+1));
I1(1:(N+1)^2,1:N*(N+1)) = kron(eye(N+1),e');

Ey = zeros((N+1)^2,N*(N+1));
for i=1:N+1
    Ey((i-1)*(N+1)+(1:N+1),:) = kron(eye(N+1),e(:,i)');
end

I1((N+1)^2+(1:(N+1)^2),N*(N+1)+(1:N*(N+1))) = Ey;
I1 = sparse(I1);


% inner-product matrix for 1-forms ( or (n-1)-forms with n=2 )
M = I1'*Qinv'*JW1*Qinv*I1;

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Inner-product matrix for 2-forms (or n-forms with n=2)                  %
% for a single spectral element                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
case 2

W2 = kron(w,w)';
JW2 = spdiags(W2./J,0,(N+1)^2,(N+1)^2);
I2 = kron(e,e)';

M = I2'*JW2*I2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end