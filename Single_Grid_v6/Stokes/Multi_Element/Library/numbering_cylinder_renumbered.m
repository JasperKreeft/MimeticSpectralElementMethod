function []=numbering_cylinder_renumbered()

% clear all
% close all
% clc
% 
% N=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

N2=N*N;

% numbering 0-cells
globalnr_0 = zeros((N+1)^2,12);

el1 = reshape((1:(N+1)^2)',N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^2,1);
ind = max(globalnr_0(:,1));

el2 = [ reshape(ind+(1:N*(N+1)),N,N+1) ; el1(1,:) ];
globalnr_0(:,2) = reshape(el2,(N+1)^2,1);
ind = max(globalnr_0(:,2));

el3 = [ el1(:,N+1) reshape(ind+(1:N*(N+1)),N+1,N) ];
globalnr_0(:,3) = reshape(el3,(N+1)^2,1);
ind = max(globalnr_0(:,3));

el4 = [ el2(:,N+1) [ reshape(ind+(1:N*N)',N,N) ; el3(1,2:N+1) ] ];
globalnr_0(:,4) = reshape(el4,(N+1)^2,1);
ind = max(globalnr_0(:,4));

el5 = [ el3(N+1,:) ; reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,5) = reshape(el5,(N+1)^2,1);
ind = max(globalnr_0(:,5));

el6 = [ el5(N+1,:) ; reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,6) = reshape(el6,(N+1)^2,1);
ind = max(globalnr_0(:,6));

el7 = [ reshape(ind+(1:N*(N+1)),N,N+1) ; el4(1,:) ];
globalnr_0(:,7) = reshape(el7,(N+1)^2,1);
ind = max(globalnr_0(:,7));

el8 = [ reshape(ind+(1:N*(N+1)),N,N+1) ; el7(1,:) ];
globalnr_0(:,8) = reshape(el8,(N+1)^2,1);
ind = max(globalnr_0(:,8));

el9 = [ el6(N+1,:) ; reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,9) = reshape(el9,(N+1)^2,1);
ind = max(globalnr_0(:,9));

el10 = [ el9(N+1,:) ; reshape(ind+(1:(N-1)*(N+1)),N-1,N+1) ; el8(1,:) ];
globalnr_0(:,10) = reshape(el10,(N+1)^2,1);
ind = max(globalnr_0(:,10));

el11 = [ reshape(ind+(1:N*(N+1)),N+1,N) el9(:,1) ];
globalnr_0(:,11) = reshape(el11,(N+1)^2,1);
ind = max(globalnr_0(:,11));

el12 = [ [ el11(N+1,1:N); reshape(ind+(1:N*N),N,N) ] el10(:,1) ];
globalnr_0(:,12) = reshape(el12,(N+1)^2,1);
ind = max(globalnr_0(:,12));


% numbering 1-cells

globalnr_1v = zeros(N*(N+1),12);
globalnr_1h = zeros(N*(N+1),12);

el1v = reshape((1:N*(N+1))',N+1,N);
globalnr_1v(:,1) = reshape(el1v,N*(N+1),1);
ind = max(globalnr_1v(:,1));

el1h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,1) = reshape(el1h',N*(N+1),1);
ind = max(globalnr_1h(:,1));

el2v = [ reshape(ind+(1:N2),N,N) ; el1v(1,:) ];
globalnr_1v(:,2) = reshape(el2v,N*(N+1),1);
ind = max(globalnr_1v(:,2));

el2h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,2) = reshape(el2h',N*(N+1),1);
ind = max(globalnr_1h(:,2));

el3v = ind + el1v;
globalnr_1v(:,3) = reshape(el3v,N*(N+1),1);
ind = max(globalnr_1v(:,3));

el3h = [ el1h(:,N+1) reshape(ind+(1:N*N),N,N)' ];
globalnr_1h(:,3) = reshape(el3h',N*(N+1),1);
ind = max(globalnr_1h(:,3));

el4v = [ reshape(ind+(1:N2),N,N) ; el3v(1,:) ];
globalnr_1v(:,4) = reshape(el4v,N*(N+1),1);
ind = max(globalnr_1v(:,4));

el4h = [ el2h(:,N+1) reshape(ind+(1:N*N)',N,N)' ];
globalnr_1h(:,4) = reshape(el4h',N*(N+1),1);
ind = max(globalnr_1h(:,4));

el5v = [ el3v(N+1,:) ; reshape(ind+(1:N*N),N,N) ];
globalnr_1v(:,5) = reshape(el5v,N*(N+1),1);
ind = max(globalnr_1v(:,5));

el5h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,5) = reshape(el5h',N*(N+1),1);
ind = max(globalnr_1h(:,5));

el6v = [ el5v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,6) = reshape(el6v,N*(N+1),1);
ind = max(globalnr_1v(:,6));

el6h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,6) = reshape(el6h',N*(N+1),1);
ind = max(globalnr_1h(:,6));

el7v = [ reshape(ind+(1:N*N),N,N) ; el4v(1,:) ];
globalnr_1v(:,7) = reshape(el7v,N*(N+1),1);
ind = max(globalnr_1v(:,7));

el7h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,7) = reshape(el7h',N*(N+1),1);
ind = max(globalnr_1h(:,7));

el8v = [ reshape(ind+(1:N2)',N,N) ; el7v(1,:) ];
globalnr_1v(:,8) = reshape(el8v,N*(N+1),1);
ind = max(globalnr_1v(:,8));

el8h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,8) = reshape(el8h',N*(N+1),1);
ind = max(globalnr_1h(:,8));

el9v = [ el6v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,9) = reshape(el9v,N*(N+1),1);
ind = max(globalnr_1v(:,9));

el9h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,9) = reshape(el9h',N*(N+1),1);
ind = max(globalnr_1h(:,9));

el10v = [ el9v(N+1,:) ; reshape(ind+(1:(N-1)*N),N,N-1)' ; el8v(1,:) ];
globalnr_1v(:,10) = reshape(el10v,N*(N+1),1);
ind = max(globalnr_1v(:,10));

el10h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,10) = reshape(el10h',N*(N+1),1);
ind = max(globalnr_1h(:,10));

el11v = reshape(ind+(1:N*(N+1)),N+1,N);
globalnr_1v(:,11) = reshape(el11v,N*(N+1),1);
ind = max(globalnr_1v(:,11));

el11h = [ reshape(ind+(1:N*N),N,N)' el9h(:,1) ];
globalnr_1h(:,11) = reshape(el11h',N*(N+1),1);
ind = max([ globalnr_1h(:,11);globalnr_1v(:,11)]);

el12v = [ el11v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,12) = reshape(el12v,N*(N+1),1);
ind = max(globalnr_1v(:,12));

el12h = [ reshape(ind+(1:N2),N,N)' el10h(:,1) ];
globalnr_1h(:,12) = reshape(el12h',N*(N+1),1);
ind = max(globalnr_1h(:,12));

% numbering 2-cells

globalnr_2 = zeros(N^2,12);
for i=1:12
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_0 = max(max(globalnr_0));
nr_1 = max(max([globalnr_1v globalnr_1h]));
nr_2 = max(max(globalnr_2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%