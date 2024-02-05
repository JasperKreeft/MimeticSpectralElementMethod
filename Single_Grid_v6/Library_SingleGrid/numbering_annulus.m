function []=numbering_annulus()

global N
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

N2=N*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_0 = zeros((N+1)^2,4);

el1 = reshape((1:(N+1)^2)',N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^2,1);
ind = max(globalnr_0(:,1));

el2 = [ (N+1):(N+1):(N+1)^2
        reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,2) = reshape(el2,(N+1)^2,1);
ind = max(globalnr_0(:,2));

el3 = [ el2(N+1,:)
        reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,3) = reshape(el3,(N+1)^2,1);
ind = max(globalnr_0(:,3));

el4 = [ el3(N+1,:)
        reshape(ind+(1:(N-1)*(N+1)),N-1,N+1)
        el1(1,:) ];
globalnr_0(:,4) = reshape(el4,(N+1)^2,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_1v = zeros(N*(N+1),4);
globalnr_1h = zeros(N*(N+1),4);

el1v = reshape((1:N*(N+1)),N+1,N);
globalnr_1v(:,1) = reshape(el1v,N*(N+1),1);
ind = max(globalnr_1v(:,1));

el1h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,1) = reshape(el1h',N*(N+1),1);
ind = max(globalnr_1h(:,1));

el2v = [ el1v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,2) = reshape(el2v,N*(N+1),1);
ind = max(globalnr_1v(:,2));

el2h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,2) = reshape(el2h',N*(N+1),1);
ind = max(globalnr_1h(:,2));

el3v = [ el2v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,3) = reshape(el3v,N*(N+1),1);
ind = max(globalnr_1v(:,3));

el3h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,3) = reshape(el3h',N*(N+1),1);
ind = max(globalnr_1h(:,3));

el4v = [ el3v(N+1,:) ; reshape(ind+(1:(N-1)*N)',N-1,N) ; el1v(1,:) ];
globalnr_1v(:,4) = reshape(el4v,N*(N+1),1);
ind = max(globalnr_1v(:,4));

el4h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,4) = reshape(el4h',N*(N+1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_2 = zeros(N^2,4);
for i=1:4
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_0 = max(max(globalnr_0));
nr_1 = max(max(globalnr_1h));
nr_2 = max(max(globalnr_2));