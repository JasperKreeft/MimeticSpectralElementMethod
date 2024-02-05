function varargout = force_Lshape(X,Y)

global N
global globalnr_1h globalnr_1v
global nr_1



f1 = zeros(4*N,N+1);
R = [1 1 2];
C = [1 2 1];
for rc = 1:3
    r = R(rc);
    c = C(rc);
    ind = globalnr_1h((1:N)+(c-1)*N,(r-1)*N+(r>1)+(1:N+(r==1)));
    f1(ind) = diff(X);
end

f2 = zeros(4*N+1,N);

F = zeros(nr_1,1);
F(globalnr_1h) = f1;
F(globalnr_1v) = f2;

if nargout==1
    varargout = {F};
elseif nargout==2
    varargout(1) = {f1};
    varargout(2) = {f2};
end