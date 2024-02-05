function []=numbering(testcase)

% clc
% testcase = 'square_periodic';
% N=1;
% numRows = 3;
% numColumns = 3;
% numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N numElements numRows numColumns
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

if isempty(numElements); numElements = 1; end
if isempty(numRows); numRows = sqrt(numElements); numColumns = numRows; end

N2=N*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch testcase
    
%% Square %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'square'
        
% numbering 0-cells
globalnr_0 = zeros((N+1)^2,numElements);

globalnr_0(:,1) = 1:(N+1)^2;

inda = 1:N+1:(N+1)^2;
indb = N+1:N+1:(N+1)^2;
indc = 1:(N+1)^2; indc(inda) = [];
for c=2:numColumns
    globalnr_0(inda,c) = globalnr_0(indb,c-1);
    globalnr_0(indc,c) = (c-1)*N*(N+1)+(N+1)+(1:N*(N+1));
end

indd = 1:N+1;
inde = N*(N+1)+1:(N+1)^2;
indg = 1:(N+1)^2; indg(inda) = []; indg(1:N) = [];
for r=2:numRows
    indf = (r-1)*numColumns*N*N + 1 + numColumns*N + (r-1)*N;
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        globalnr_0(indd,rc) = globalnr_0(inde,rc-numColumns);
        if c==1
            globalnr_0(N+2:(N+1)^2,rc) = indf + (1:N*(N+1));
            indf = indf + N*(N+1);
        else
            globalnr_0(inda,rc) = globalnr_0(indb,rc-1);
            globalnr_0(indg,rc) = indf+(1:N^2);
            indf = indf + N^2;
        end
    end
end


% numbering 1-cells
globalnr_1h = zeros(N*(N+1),numElements);
globalnr_1v = zeros(N*(N+1),numElements);

ind = 0;
inda = 1:N+1:N*(N+1);
indb = N+1:N+1:N*(N+1);
indc = 1:N*(N+1); indc(inda) = [];
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        if c==1
globalnr_1v(:,rc) = ind + (1:N*(N+1));
ind = ind + N*(N+1);
        else
globalnr_1v(inda,rc) = globalnr_1v(indb,rc-1);
globalnr_1v(indc,rc) = ind+(1:N2);
ind = ind + N2;
        end
        if r==1
globalnr_1h(:,rc) = ind+(1:N*(N+1));
ind = ind + N*(N+1);            
        else
globalnr_1h(inda,rc) = globalnr_1h(indb,rc-numColumns);
globalnr_1h(indc,rc) = ind+(1:N2);
ind = ind + N2;
        end
    end
end


% numbering 2-cells
globalnr_2 = zeros(N^2,numElements);
for i=1:numElements
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end

%% Square periodic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'square_periodic'

% numbering 0-cells

if numElements==1
    
    globalnr_0 = zeros((N+1)^2,1);
    for i=1:N
        globalnr_0((i-1)*(N+1)+(1:N),1) = (i-1)*N+(1:N);
        globalnr_0(N+1:N+1:(N+1)^2,1) = globalnr_0(1:N+1:(N+1)^2,1);
        globalnr_0(N*(N+1)+1:(N+1)^2,1) = globalnr_0(1:N+1,1);
    end
    
else
    
globalnr_0 = zeros((N+1)^2,numElements);
globalnr_0(:,1) = 1:(N+1)^2;

indi = 1:(N+1)^2; indi([1:N+1:(N+1)^2 N+1:N+1:(N+1)^2]) = [];
indj = 1:(N+1)^2; indj([1:N+1 1:N+1:(N+1)^2 N+1:N+1:(N+1)^2]) = [];
inda = 1:N+1:(N+1)^2;
indb = N+1:N+1:(N+1)^2;
indc = 1:(N+1)^2; indc(inda) = [];
for c=2:numColumns-1
    globalnr_0(inda,c) = globalnr_0(indb,c-1);
    globalnr_0(indc,c) = (c-1)*N*(N+1)+(N+1)+(1:N*(N+1));
end
globalnr_0(inda,numColumns) = globalnr_0(indb,numColumns-1);
globalnr_0(indi,numColumns) = (numColumns-1)*N*(N+1) +(N+1) + (1:(N-1)*(N+1));
globalnr_0(indb,numColumns) = globalnr_0(inda,1);

indd = 1:N+1;
inde = N*(N+1)+1:(N+1)^2;
indg = 1:(N+1)^2; indg(inda) = []; indg(1:N) = [];
indk = 1:(N+1)^2; indk([inda indd inde]) = [];
indl = 1:(N+1)^2; indl([inda indb indd inde]) = [];
indf = numColumns*N*(N+1);
for r=2:numRows-1
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        globalnr_0(indd,rc) = globalnr_0(inde,rc-numColumns);
        if c==1
            globalnr_0(N+2:(N+1)^2,rc) = indf + (1:N*(N+1));
            indf = indf + N*(N+1);
        elseif c==numColumns
            globalnr_0(indd,r*numColumns) = globalnr_0(inde,(r-1)*numColumns);
        globalnr_0(inda,r*numColumns) = globalnr_0(indb,r*numColumns-1);
        globalnr_0(indj,r*numColumns) = indf + (1:(N-1)*N);
        globalnr_0(indb,r*numColumns) = globalnr_0(inda,(r-1)*numColumns+1);
        indf = indf + (N-1)*N;
        else
            globalnr_0(inda,rc) = globalnr_0(indb,rc-1);
            globalnr_0(indg,rc) = indf+(1:N^2);
            indf = indf + N^2;
        end
    end
end
for c=1:numColumns
    rc = c+(numRows-1)*numColumns;
    globalnr_0(indd,rc) = globalnr_0(inde,rc-numColumns);
    if c==1
        globalnr_0(N+2:N*(N+1),rc) = indf + (1:(N-1)*(N+1));
        indf = indf + (N-1)*(N+1);
    elseif c==numColumns
        globalnr_0(indd,numElements) = globalnr_0(inde,(numRows-1)*numColumns);
        globalnr_0(inda,numElements) = globalnr_0(indb,numElements-1);
        globalnr_0(inde,numElements) = globalnr_0(indd,numColumns);
        globalnr_0(indb,numElements) = globalnr_0(inda,(numRows-1)*numColumns+1);
        globalnr_0(indl,numElements) = indf + (1:(N-1)^2);
    else
        globalnr_0(inda,rc) = globalnr_0(indb,rc-1);
        globalnr_0(indk,rc) = indf+(1:(N-1)*N);
        indf = indf + (N-1)*N;
    end
    globalnr_0(inde,rc) = globalnr_0(indd,c);
end

end

% numbering 1-cells
globalnr_1h = zeros(N*(N+1),numElements);
globalnr_1v = zeros(N*(N+1),numElements);


inda = 1:N+1:N*(N+1);
indb = N+1:N+1:N*(N+1);
indc = 1:N*(N+1); indc(inda) = [];
indd = 1:N*(N+1); indd([inda indb]) = [];
inde = 1:N*(N+1); inde(indb) = [];

if numElements==1

globalnr_1v(inde,1) = 1:N2;
globalnr_1v(indb,1) = globalnr_1v(inda,1);
globalnr_1h(inde,1) = N2+(1:N2);
globalnr_1h(indb,1) = globalnr_1h(inda,1);

else

ind = 0;
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        if c==1
globalnr_1v(:,rc) = ind + (1:N*(N+1));
ind = ind + N*(N+1);
        elseif c==numColumns
globalnr_1v(inda,rc) = globalnr_1v(indb,rc-1);
globalnr_1v(indd,rc) = ind+(1:(N-1)*N);
globalnr_1v(indb,rc) = globalnr_1v(inda,rc-numColumns+1);
ind = ind + (N-1)*N;
        else
globalnr_1v(inda,rc) = globalnr_1v(indb,rc-1);
globalnr_1v(indc,rc) = ind+(1:N2);
ind = ind + N2;
        end
        if r==1
globalnr_1h(:,rc) = ind+(1:N*(N+1));
ind = ind + N*(N+1);            
        elseif r==numRows
globalnr_1h(inda,rc) = globalnr_1h(indb,rc-numColumns);
globalnr_1h(indd,rc) = ind+(1:N*(N-1));
globalnr_1h(indb,rc) = globalnr_1h(inda,c);
ind = ind + N*(N-1);
        else
globalnr_1h(inda,rc) = globalnr_1h(indb,rc-numColumns);
globalnr_1h(indc,rc) = ind+(1:N2);
ind = ind + N2;
        end
    end
end

end


% numbering 2-cells
globalnr_2 = zeros(N^2,numElements);
for i=1:numElements
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end



%% Triangle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'triangle'
        
% numbering 0-cells
element1 = (1:(N+1)^2)';
el2 = [ (N+1):(N+1):(N+1)^2
  reshape((N+1)^2+(1:N*(N+1)),N,N+1) ];
element2 = reshape(el2,(N+1)^2,1);

el3 = [ (N+1)^2-(0:N)
        [ (N+1)^2+N^2+(1:N)' reshape((N+1)^2+N*(N+1)+(1:N*N)',N,N)] ];
element3  = reshape(el3,(N+1)^2,1);
         
globalnr_0 = [ element1 element2 element3 ];

% numbering 1-cells
element1 = (1:N*(N+1))';
el2 = [ N+1:N+1:N*(N+1)  ;  reshape(2*N*(N+1)+(1:N^2),N,N) ];
element2 = reshape(el2,N*(N+1),1);
el3 = [ 2*N*(N+1):-(N+1):N*(N+1)+1 ; reshape(3*N*(N+1)+N^2+(1:N^2),N,N) ];
element3 = reshape(el3,N*(N+1),1);

globalnr_1v = [ element1 element2 element3 ];

element1 = (N*(N+1)+(1:N*(N+1)))';
element2 = 2*N*(N+1)+N^2+(1:N*(N+1))';
el3 = [ 2*N*(N+1)+N^2+(N+1:N+1:N*(N+1)) ; reshape(3*N*(N+1)+2*N^2+(1:N^2),N,N) ]';
element3 = reshape(el3',N*(N+1),1);

globalnr_1h = [ element1 element2 element3 ];
            
% numbering 2-cells
globalnr_2 = zeros(N*N,3);
for i=1:3
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end

%% Triangle E9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'triangle_E9'

% numbering 0-cells
globalnr_0 = zeros((N+1)^2,9);

el1 = reshape((1:(N+1)^2)',N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^2,1);
ind = max(globalnr_0(:,1));

el2 = [ el1(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,2) = reshape(el2,(N+1)^2,1);
ind = max(globalnr_0(:,2));

el3 = [ el2(N+1,:)
       reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,3) = reshape(el3,(N+1)^2,1);
ind = max(globalnr_0(:,3));

el4 = [ el3(N+1,:)
       reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,4) = reshape(el4,(N+1)^2,1);
ind = max(globalnr_0(:,4));

el5 = [ el1(N+1:-1:1,N+1)'
        el2(2:N+1,N+1) reshape(ind+(1:N*N),N,N) ];
globalnr_0(:,5) = reshape(el5,(N+1)^2,1);
ind = max(globalnr_0(:,5));

el6 = [ el4(:,N+1) [ el3(N:-1:1,N+1)' ; reshape(ind+(1:N*N),N,N) ] ];
globalnr_0(:,6) = reshape(el6,(N+1)^2,1);
ind = max(globalnr_0(:,6));

el7 = [ el5(N+1,N+1:-1:1)' reshape(ind+(1:N*(N+1)),N+1,N) ];
globalnr_0(:,7) = reshape(el7,(N+1)^2,1);
ind = max(globalnr_0(:,7));

el8 = [ el6(:,N+1) [ el7(N+1,2:N+1) ; reshape(ind+(1:N*N),N,N) ] ];
globalnr_0(:,8) = reshape(el8,(N+1)^2,1);
ind = max(globalnr_0(:,8));

el9 = [ el8(:,N+1) [ el7(N:-1:1,N+1)' ; reshape(ind+(1:N*N),N,N) ] ];
globalnr_0(:,9) = reshape(el9,(N+1)^2,1);


% numbering 1-cells
globalnr_1v = zeros(N*(N+1),9);
globalnr_1h = zeros(N*(N+1),9);

el1v = reshape((1:N*(N+1))',N+1,N);
globalnr_1v(:,1) = reshape(el1v,N*(N+1),1);
ind = max(globalnr_1v(:,1));

el1h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,1) = reshape(el1h',N*(N+1),1);
ind = max(globalnr_1h(:,1));

el2v = [ el1v(N+1,:) ; reshape(ind+(1:N*N)',N,N) ];
globalnr_1v(:,2) = reshape(el2v,N*(N+1),1);
ind = max(globalnr_1v(:,2));

el2h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,2) = reshape(el2h',N*(N+1),1);
ind = max(globalnr_1h(:,2));

el3v = [ el2v(N+1,:) ; reshape(ind+(1:N*N)',N,N) ];
globalnr_1v(:,3) = reshape(el3v,N*(N+1),1);
ind = max(globalnr_1v(:,3));

el3h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,3) = reshape(el3h',N*(N+1),1);
ind = max(globalnr_1h(:,3));

el4v = [ el3v(N+1,:) ; reshape(ind+(1:N*N)',N,N) ];
globalnr_1v(:,4) = reshape(el4v,N*(N+1),1);
ind = max(globalnr_1v(:,4));

el4h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,4) = reshape(el4h',N*(N+1),1);
ind = max(globalnr_1h(:,4));

el5v = [ el1h(N:-1:1,N+1)' ; reshape(ind+(1:N*N),N,N) ];
globalnr_1v(:,5) = reshape(el5v,N*(N+1),1);
ind = max(globalnr_1v(:,5));

el5h = [ el2h(:,N+1) reshape(ind+(1:N*N),N,N)' ];
globalnr_1h(:,5) = reshape(el5h',N*(N+1),1);
ind = max(globalnr_1h(:,5));

el6v = [ el3h(N:-1:1,N+1)' ; reshape(ind+(1:N*N),N,N) ];
globalnr_1v(:,6) = reshape(el6v,N*(N+1),1);
ind = max(globalnr_1v(:,6));

el6h = [ el4h(:,N+1) reshape(ind+(1:N*N),N,N)' ];
globalnr_1h(:,6) = reshape(el6h',N*(N+1),1);
ind = max(globalnr_1h(:,6));

el7v = reshape(ind+(1:N*(N+1)),N+1,N);
globalnr_1v(:,7) = reshape(el7v,N*(N+1),1);
ind = max(globalnr_1v(:,7));

el7h = [ el5v(N+1,N:-1:1)' reshape(ind+(1:N*N),N,N)' ];
globalnr_1h(:,7) = reshape(el7h',N*(N+1),1);
ind = max(globalnr_1h(:,7));

el8v = [ el7v(N+1,:) ; reshape(ind+(1:N*N),N,N) ];
globalnr_1v(:,8) = reshape(el8v,N*(N+1),1);
ind = max(globalnr_1v(:,8));

el8h = [ el6h(:,N+1) reshape(ind+(1:N*N),N,N)' ];
globalnr_1h(:,8) = reshape(el8h',N*(N+1),1);
ind = max(globalnr_1h(:,8));


el9v = [ el7h(N:-1:1,N+1)' ; reshape(ind+(1:N*N),N,N) ];
globalnr_1v(:,9) = reshape(el9v,N*(N+1),1);
ind = max(globalnr_1v(:,9));

el9h = [ el8h(:,N+1) reshape(ind+(1:N*N),N,N)' ];
globalnr_1h(:,9) = reshape(el9h',N*(N+1),1);

% numbering 2-cells

globalnr_2 = zeros(N*N,9);

for i=1:9
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end

%% Semicylinder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'semicylinder'

% numbering 0-cells
globalnr_0 = zeros((N+1)^2,6);

el1 = reshape((1:(N+1)^2)',N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^2,1);
ind = max(globalnr_0(:,1));

el2 = [ (N+1):(N+1):(N+1)^2
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,2) = reshape(el2,(N+1)^2,1);
ind = max(globalnr_0(:,2));

el3 = [ el2(:,1)' ; reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,3) = reshape(el3,(N+1)^2,1);
ind = max(globalnr_0(:,3));

el4 = [ el3(N+1,:) ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,4) = reshape(el4,(N+1)^2,1);
ind = max(globalnr_0(:,4));

el5 = [ fliplr(el4(N+1,:))' reshape(ind+(1:N*(N+1))',N+1,N) ];
globalnr_0(:,5) = reshape(el5,(N+1)^2,1);
ind = max(globalnr_0(:,5));

el6 = [ el5(N+1,:) ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,6) = reshape(el6,(N+1)^2,1);


% numbering 1-cells
globalnr_1v = zeros(N*(N+1),6);
globalnr_1h = zeros(N*(N+1),6);

el1v = reshape((1:N*(N+1))',N+1,N);
globalnr_1v(:,1) = reshape(el1v,N*(N+1),1);
ind = max(globalnr_1v(:,1));

el1h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,1) = reshape(el1h',N*(N+1),1);
ind = max(globalnr_1h(:,1));

el2v = [ el1v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,2) = reshape(el2v,N*(N+1),1);
ind = max(globalnr_1v(:,2));

el2h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,2) = reshape(el2h',N*(N+1),1);
ind = max(globalnr_1h(:,2));

el3v = [ el2h(:,1)' ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,3) = reshape(el3v,N*(N+1),1);
ind = max(globalnr_1v(:,3));

el3h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,3) = reshape(el3h',N*(N+1),1);
ind = max(globalnr_1h(:,3));

el4v = [ el3v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,4) = reshape(el4v,N*(N+1),1);
ind = max(globalnr_1v(:,4));

el4h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,4) = reshape(el4h',N*(N+1),1);
ind = max(globalnr_1h(:,4));

el5v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,5) = reshape(el5v,N*(N+1),1);
ind = max(globalnr_1v(:,5));

el5h = [ fliplr(el4v(N+1,:)) ; reshape(ind+(1:N2)',N,N)  ]';
globalnr_1h(:,5) = reshape(el5h',N*(N+1),1);
ind = max(globalnr_1h(:,5));

el6v = [ el5v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,6) = reshape(el6v,N*(N+1),1);
ind = max(globalnr_1v(:,6));

el6h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,6) = reshape(el6h',N*(N+1),1);

% numbering 2-cells
globalnr_2 = zeros(N^2,6);
for i=1:6
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end



    case 'semicylinder_v2'
        

% numbering 0-cells
globalnr_0 = zeros((N+1)^2,6);

el1 = reshape((1:(N+1)^2)',N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^2,1);
ind = max(globalnr_0(:,1));

el2 = [ N*(N+1)+1:(N+1)^2
        reshape(ind+(1:N*(N+1)),N+1,N)' ]';
globalnr_0(:,2) = reshape(el2,(N+1)^2,1);
ind = max(globalnr_0(:,2));

el3 = [ el2(N+1,:) ; reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,3) = reshape(el3,(N+1)^2,1);
ind = max(globalnr_0(:,3));

el4 = [ el3(N+1,:) ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,4) = reshape(el4,(N+1)^2,1);
ind = max(globalnr_0(:,4));

el5 = [ el4(N+1,:) ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,5) = reshape(el5,(N+1)^2,1);
ind = max(globalnr_0(:,5));

el6 = [ reshape(ind+(1:N*(N+1)),N+1,N) el5(:,1) ];
globalnr_0(:,6) = reshape(el6,(N+1)^2,1);


% numbering 1-cells
globalnr_1v = zeros(N*(N+1),6);
globalnr_1h = zeros(N*(N+1),6);

el1v = reshape((1:N*(N+1)),N+1,N);
globalnr_1v(:,1) = reshape(el1v,N*(N+1),1);
ind = max(globalnr_1v(:,1));

el1h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,1) = reshape(el1h',N*(N+1),1);
ind = max(globalnr_1h(:,1));

el2v = reshape(ind+(1:N*(N+1)),N+1,N);
globalnr_1v(:,2) = reshape(el2v,N*(N+1),1);
ind = max(globalnr_1v(:,2));

el2h = [el1h(:,N+1) reshape(ind+(1:N2),N,N)'];
globalnr_1h(:,2) = reshape(el2h',N*(N+1),1);
ind = max(globalnr_1h(:,2));

el3v = [ el2v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,3) = reshape(el3v,N*(N+1),1);
ind = max(globalnr_1v(:,3));

el3h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,3) = reshape(el3h',N*(N+1),1);
ind = max(globalnr_1h(:,3));

el4v = [ el3v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,4) = reshape(el4v,N*(N+1),1);
ind = max(globalnr_1v(:,4));

el4h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,4) = reshape(el4h',N*(N+1),1);
ind = max(globalnr_1h(:,4));

el5v = [ el4v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,5) = reshape(el5v,N*(N+1),1);
ind = max(globalnr_1v(:,5));

el5h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,5) = reshape(el5h',N*(N+1),1);
ind = max(globalnr_1h(:,5));

el6v = reshape(ind+(1:N*(N+1)),N+1,N);
globalnr_1v(:,6) = reshape(el6v,N*(N+1),1);
ind = max(globalnr_1v(:,6));

el6h = [ reshape(ind+(1:N2),N,N)' el5h(:,1) ];
globalnr_1h(:,6) = reshape(el6h',N*(N+1),1);

% numbering 2-cells
globalnr_2 = zeros(N^2,6);
for i=1:6
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end


%% Cylinder %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'cylinder'

% numbering 0-cells
globalnr_0 = zeros((N+1)^2,12);

el1 = reshape((1:(N+1)^2)',N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^2,1);
ind = max(globalnr_0(:,1));

el2 = [ (N+1):(N+1):(N+1)^2
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,2) = reshape(el2,(N+1)^2,1);
ind = max(globalnr_0(:,2));

el3 = [ el2(:,1)' ; reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,3) = reshape(el3,(N+1)^2,1);
ind = max(globalnr_0(:,3));

el4 = [ el3(N+1,:) ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,4) = reshape(el4,(N+1)^2,1);
ind = max(globalnr_0(:,4));

el5 = [ fliplr(el4(N+1,:))' reshape(ind+(1:N*(N+1))',N+1,N) ];
globalnr_0(:,5) = reshape(el5,(N+1)^2,1);
ind = max(globalnr_0(:,5));

el6 = [ el5(N+1,:) ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,6) = reshape(el6,(N+1)^2,1);
ind = max(globalnr_0(:,6));

el7 = [ el1(:,N+1) reshape(ind+(1:N*(N+1))',N+1,N) ];
globalnr_0(:,7) = reshape(el7,(N+1)^2,1);
ind = max(globalnr_0(:,7));

el8 = [ el2(:,N+1) [ el7(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ] ];
globalnr_0(:,8) = reshape(el8,(N+1)^2,1);
ind = max(globalnr_0(:,8));

el9 = [ fliplr(el8(:,N+1)') ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,9) = reshape(el9,(N+1)^2,1);
ind = max(globalnr_0(:,9));

el10 = [ el9(N+1,:) ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,10) = reshape(el10,(N+1)^2,1);
ind = max(globalnr_0(:,10));

el11 = [ el5(:,N+1) reshape(ind+(1:(N-1)*(N+1))',N+1,N-1) el10(N+1,:)' ];
globalnr_0(:,11) = reshape(el11,(N+1)^2,1);
ind = max(globalnr_0(:,11));

el12 = [ el6(:,N+1) [ el11(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ] ];
globalnr_0(:,12) = reshape(el12,(N+1)^2,1);
ind = max(globalnr_0(:,12));

% numbering 1-cells

globalnr_1v = zeros(N*(N+1),12);
globalnr_1h = zeros(N*(N+1),12);

el1v = reshape((1:N*(N+1))',N+1,N);
globalnr_1v(:,1) = reshape(el1v,N*(N+1),1);
ind = max(globalnr_1v(:,1));

el1h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,1) = reshape(el1h',N*(N+1),1);
ind = max(globalnr_1h(:,1));

el2v = [ el1v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,2) = reshape(el2v,N*(N+1),1);
ind = max(globalnr_1v(:,2));

el2h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,2) = reshape(el2h',N*(N+1),1);
ind = max(globalnr_1h(:,2));

el3v = [ el2h(:,1)' ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,3) = reshape(el3v,N*(N+1),1);
ind = max(globalnr_1v(:,3));

el3h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,3) = reshape(el3h',N*(N+1),1);
ind = max(globalnr_1h(:,3));

el4v = [ el3v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,4) = reshape(el4v,N*(N+1),1);
ind = max(globalnr_1v(:,4));

el4h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,4) = reshape(el4h',N*(N+1),1);
ind = max(globalnr_1h(:,4));

el5v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,5) = reshape(el5v,N*(N+1),1);
ind = max(globalnr_1v(:,5));

el5h = [ fliplr(el4v(N+1,:)) ; reshape(ind+(1:N2)',N,N)  ]';
globalnr_1h(:,5) = reshape(el5h',N*(N+1),1);
ind = max(globalnr_1h(:,5));

el6v = [ el5v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,6) = reshape(el6v,N*(N+1),1);
ind = max(globalnr_1v(:,6));

el6h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,6) = reshape(el6h',N*(N+1),1);
ind = max(globalnr_1h(:,6));

el7v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,7) = reshape(el7v,N*(N+1),1);
ind = max(globalnr_1v(:,7));

el7h = [ el1h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,7) = reshape(el7h',N*(N+1),1);
ind = max(globalnr_1h(:,7));

el8v = [ el7v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,8) = reshape(el8v,N*(N+1),1);
ind = max(globalnr_1v(:,8));

el8h = [ el2h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,8) = reshape(el8h',N*(N+1),1);
ind = max(globalnr_1h(:,8));

el9v = [ fliplr(el8h(:,N+1)') ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,9) = reshape(el9v,N*(N+1),1);
ind = max(globalnr_1v(:,9));

el9h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,9) = reshape(el9h',N*(N+1),1);
ind = max(globalnr_1h(:,9));

el10v = [ el9v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,10) = reshape(el10v,N*(N+1),1);
ind = max(globalnr_1v(:,10));

el10h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,10) = reshape(el10h',N*(N+1),1);
ind = max(globalnr_1h(:,10));

el11v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,11) = reshape(el11v,N*(N+1),1);
ind = max(globalnr_1v(:,11));

el11h = [ el5h(:,N+1) reshape(ind+(1:(N-1)*N)',N-1,N)' el10v(N+1,:)' ];
globalnr_1h(:,11) = reshape(el11h',N*(N+1),1);
ind = max([ globalnr_1h(:,11);globalnr_1v(:,11)]);

el12v = [ el11v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,12) = reshape(el12v,N*(N+1),1);
ind = max(globalnr_1v(:,12));

el12h = [ el6h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,12) = reshape(el12h',N*(N+1),1);
ind = max(globalnr_1h(:,12));

% numbering 2-cells

globalnr_2 = zeros(N^2,12);
for i=1:12
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end

%% Cylinder E20 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'cylinder_E20'

% 0-cells
globalnr_0 = zeros((N+1)^2,20);

el1 = reshape((1:(N+1)^2)',N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^2,1);
ind = max(globalnr_0(:,1));

el2 = [ el1(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,2) = reshape(el2,(N+1)^2,1);
ind = max(globalnr_0(:,2));

el3 = [ el2(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,3) = reshape(el3,(N+1)^2,1);
ind = max(globalnr_0(:,3));

el4 = [ el3(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,4) = reshape(el4,(N+1)^2,1);
ind = max(globalnr_0(:,4));

el5 = [ el1(:,N+1) ind+reshape(1:N*(N+1),N+1,N) ];
globalnr_0(:,5) = reshape(el5,(N+1)^2,1);
ind = max(globalnr_0(:,5));

el6 = [ el5(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,6) = reshape(el6,(N+1)^2,1);
ind = max(globalnr_0(:,6));

el7 = [ el2(:,N+1) [ el6(2:N+1,1)' ; reshape(ind+(1:N2),N,N) ]];
globalnr_0(:,7) = reshape(el7,(N+1)^2,1);
ind = max(globalnr_0(:,7));

el8 = [ el3(:,N+1) [ el7(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ]];
globalnr_0(:,8) = reshape(el8,(N+1)^2,1);
ind = max(globalnr_0(:,8));

el9 = [ fliplr(el8(N+1,:))' reshape(ind+(1:N*(N+1))',N+1,N) ];
globalnr_0(:,9) = reshape(el9,(N+1)^2,1);
ind = max(globalnr_0(:,9));

el10 = [ el4(:,N+1) [ el9(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ]];
globalnr_0(:,10) = reshape(el10,(N+1)^2,1);
ind = max(globalnr_0(:,10));

el11 = [ el5(:,N+1) reshape(ind+(1:N*(N+1))',N+1,N) ];
globalnr_0(:,11) = reshape(el11,(N+1)^2,1);
ind = max(globalnr_0(:,11));

el12 = [ el6(:,N+1) [ el11(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ] ];
globalnr_0(:,12) = reshape(el12,(N+1)^2,1);
ind = max(globalnr_0(:,12));

el13 = [ fliplr(el12(:,N+1)') ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,13) = reshape(el13,(N+1)^2,1);
ind = max(globalnr_0(:,13));

el14 = [ el13(N+1,:) ; reshape(ind+(1:N*(N+1))',N,N+1) ];
globalnr_0(:,14) = reshape(el14,(N+1)^2,1);
ind = max(globalnr_0(:,14));

el15 = [ el9(:,N+1) reshape(ind+(1:(N-1)*(N+1))',N+1,N-1) el14(N+1,:)' ];
globalnr_0(:,15) = reshape(el15,(N+1)^2,1);
ind = max(globalnr_0(:,15));

el16 = [ el10(:,N+1) [ el15(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ] ];
globalnr_0(:,16) = reshape(el16,(N+1)^2,1);
ind = max(globalnr_0(:,16));

el17 = [ el11(:,N+1) reshape(ind+(1:N*(N+1)),N+1,N) ];
globalnr_0(:,17) = reshape(el17,(N+1)^2,1);
ind = max(globalnr_0(:,17));


el18 = [ el13(:,N+1) [ el17(N+1,2:N+1) ; reshape(ind+(1:N2),N,N) ]];
globalnr_0(:,18) = reshape(el18,(N+1)^2,1);
ind = max(globalnr_0(:,18));

el19 = [ el14(:,N+1) [ el18(N+1,2:N+1) ; reshape(ind+(1:N2),N,N) ]];
globalnr_0(:,19) = reshape(el19,(N+1)^2,1);
ind = max(globalnr_0(:,19));

el20 = [ el16(:,N+1) [ el19(N+1,2:N+1) ; reshape(ind+(1:N2),N,N) ]];
globalnr_0(:,20) = reshape(el20,(N+1)^2,1);


% 1-cells
globalnr_1v = zeros(N*(N+1),20);
globalnr_1h = zeros(N*(N+1),20);

el1v = reshape((1:N*(N+1))',N+1,N);
globalnr_1v(:,1) = reshape(el1v,N*(N+1),1);
ind = max(globalnr_1v(:,1));

el1h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,1) = reshape(el1h',N*(N+1),1);
ind = max(globalnr_1h(:,1));

el2v = [ el1v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,2) = reshape(el2v,N*(N+1),1);
ind = max(globalnr_1v(:,2));

el2h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,2) = reshape(el2h',N*(N+1),1);
ind = max(globalnr_1h(:,2));

el3v = [ el2v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,3) = reshape(el3v,N*(N+1),1);
ind = max(globalnr_1v(:,3));

el3h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,3) = reshape(el3h',N*(N+1),1);
ind = max(globalnr_1h(:,3));

el4v = [ el3v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,4) = reshape(el4v,N*(N+1),1);
ind = max(globalnr_1v(:,4));

el4h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,4) = reshape(el4h',N*(N+1),1);
ind = max(globalnr_1h(:,4));

el5v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,5) = reshape(el5v,N*(N+1),1);
ind = max(globalnr_1v(:,5));

el5h = [ el1h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,5) = reshape(el5h',N*(N+1),1);
ind = max(globalnr_1h(:,5));

el6v = [ el5v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,6) = reshape(el6v,N*(N+1),1);
ind = max(globalnr_1v(:,6));

el6h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,6) = reshape(el6h',N*(N+1),1);
ind = max(globalnr_1h(:,6));

el7v = [ el6h(:,1)' ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,7) = reshape(el7v,N*(N+1),1);
ind = max(globalnr_1v(:,7));

el7h = [ el2h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,7) = reshape(el7h',N*(N+1),1);
ind = max(globalnr_1h(:,7));

el8v = [ el7v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,8) = reshape(el8v,N*(N+1),1);
ind = max(globalnr_1v(:,8));

el8h = [ el3h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,8) = reshape(el8h',N*(N+1),1);
ind = max(globalnr_1h(:,8));

el9v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,9) = reshape(el9v,N*(N+1),1);
ind = max(globalnr_1v(:,9));

el9h = [ fliplr(el8v(N+1,:)) ; reshape(ind+(1:N2)',N,N)  ]';
globalnr_1h(:,9) = reshape(el9h',N*(N+1),1);
ind = max(globalnr_1h(:,9));

el10v = [ el9v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,10) = reshape(el10v,N*(N+1),1);
ind = max(globalnr_1v(:,10));

el10h = [ el4h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,10) = reshape(el10h',N*(N+1),1);
ind = max(globalnr_1h(:,10));

el11v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,11) = reshape(el11v,N*(N+1),1);
ind = max(globalnr_1v(:,11));

el11h = [ el5h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,11) = reshape(el11h',N*(N+1),1);
ind = max(globalnr_1h(:,11));

el12v = [ el11v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,12) = reshape(el12v,N*(N+1),1);
ind = max(globalnr_1v(:,12));

el12h = [ el6h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,12) = reshape(el12h',N*(N+1),1);
ind = max(globalnr_1h(:,12));

el13v = [ fliplr(el12h(:,N+1)') ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,13) = reshape(el13v,N*(N+1),1);
ind = max(globalnr_1v(:,13));

el13h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,13) = reshape(el13h',N*(N+1),1);
ind = max(globalnr_1h(:,13));

el14v = [ el13v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,14) = reshape(el14v,N*(N+1),1);
ind = max(globalnr_1v(:,14));

el14h = reshape(ind+(1:N*(N+1))',N+1,N)';
globalnr_1h(:,14) = reshape(el14h',N*(N+1),1);
ind = max(globalnr_1h(:,14));

el15v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,15) = reshape(el15v,N*(N+1),1);
ind = max(globalnr_1v(:,15));

el15h = [ el9h(:,N+1) reshape(ind+(1:(N-1)*N)',N-1,N)' el14v(N+1,:)' ];
globalnr_1h(:,15) = reshape(el15h',N*(N+1),1);
ind = max([ globalnr_1h(:,15);globalnr_1v(:,15)]);

el16v = [ el15v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,16) = reshape(el16v,N*(N+1),1);
ind = max(globalnr_1v(:,16));

el16h = [ el10h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,16) = reshape(el16h',N*(N+1),1);
ind = max(globalnr_1h(:,16));

el17v = reshape(ind+(1:N*(N+1))',N+1,N);
globalnr_1v(:,17) = reshape(el17v,N*(N+1),1);
ind = max(globalnr_1v(:,17));

el17h = [el11h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,17) = reshape(el17h',N*(N+1),1);
ind = max(globalnr_1h(:,17));

el18v = [ el17v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,18) = reshape(el18v,N*(N+1),1);
ind = max(globalnr_1v(:,18));

el18h = [ el13h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,18) = reshape(el18h',N*(N+1),1);
ind = max(globalnr_1h(:,18));

el19v = [ el18v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,19) = reshape(el19v,N*(N+1),1);
ind = max(globalnr_1v(:,19));

el19h = [ el14h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,19) = reshape(el19h',N*(N+1),1);
ind = max(globalnr_1h(:,19));

el20v = [ el19v(N+1,:) ; reshape(ind+(1:N2)',N,N) ];
globalnr_1v(:,20) = reshape(el20v,N*(N+1),1);
ind = max(globalnr_1v(:,20));

el20h = [ el16h(:,N+1) reshape(ind+(1:N2)',N,N)' ];
globalnr_1h(:,20) = reshape(el20h',N*(N+1),1);


% 2-cells
globalnr_2 = zeros(N^2,20);
for i=1:20
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end


%% Cylinder E16 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'cylinder_E16'

% 0-cells
globalnr_0 = zeros((N+1)^2,16);

el1 = reshape((1:(N+1)^2)',N+1,N+1);
globalnr_0(:,1) = reshape(el1,(N+1)^2,1);
ind = max(globalnr_0(:,1));

el2 = [ el1(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,2) = reshape(el2,(N+1)^2,1);
ind = max(globalnr_0(:,2));

el3 = [ el2(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,3) = reshape(el3,(N+1)^2,1);
ind = max(globalnr_0(:,3));

el4 = [ el3(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,4) = reshape(el4,(N+1)^2,1);
ind = max(globalnr_0(:,4));

el5 = [ el4(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,5) = reshape(el5,(N+1)^2,1);
ind = max(globalnr_0(:,5));

el6 = [ el5(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,6) = reshape(el6,(N+1)^2,1);
ind = max(globalnr_0(:,6));

el7 = [ el6(N+1,:)
  reshape(ind+(1:N*(N+1)),N,N+1) ];
globalnr_0(:,7) = reshape(el7,(N+1)^2,1);
ind = max(globalnr_0(:,7));

el8 = [ el7(N+1,:)
  reshape(ind+(1:(N-1)*(N+1)),N-1,N+1)
         el1(1,:)     ];
globalnr_0(:,8) = reshape(el8,(N+1)^2,1);
ind = max(globalnr_0(:,8));

el9 = [ el1(:,N+1) reshape(ind+(1:N*(N+1)),N+1,N) ];
globalnr_0(:,9) = reshape(el9,(N+1)^2,1);
ind = max(globalnr_0(:,9));

el10 = [ el2(:,N+1) [ el9(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ]];
globalnr_0(:,10) = reshape(el10,(N+1)^2,1);
ind = max(globalnr_0(:,10));

el11 = [ el3(:,N+1) [ el10(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ]];
globalnr_0(:,11) = reshape(el11,(N+1)^2,1);
ind = max(globalnr_0(:,11));

el12 = [ el4(:,N+1) [ el11(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ]];
globalnr_0(:,12) = reshape(el12,(N+1)^2,1);
ind = max(globalnr_0(:,12));

el13 = [ el5(:,N+1) [ el12(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ]];
globalnr_0(:,13) = reshape(el13,(N+1)^2,1);
ind = max(globalnr_0(:,13));

el14 = [ el6(:,N+1) [ el13(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ]];
globalnr_0(:,14) = reshape(el14,(N+1)^2,1);
ind = max(globalnr_0(:,14));

el15 = [ el7(:,N+1) [ el14(N+1,2:N+1) ; reshape(ind+(1:N2)',N,N) ]];
globalnr_0(:,15) = reshape(el15,(N+1)^2,1);
ind = max(globalnr_0(:,15));

el16 = [ el8(:,N+1) [ el15(N+1,2:N+1) ; reshape(ind+(1:N*(N-1)),N-1,N) ; el9(1,2:N+1) ]];
globalnr_0(:,16) = reshape(el16,(N+1)^2,1);
ind = max(globalnr_0(:,16));

% 1-cells
globalnr_1v = zeros(N*(N+1),16);
globalnr_1h = zeros(N*(N+1),16);

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

el4v = [ el3v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,4) = reshape(el4v,N*(N+1),1);
ind = max(globalnr_1v(:,4));

el4h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,4) = reshape(el4h',N*(N+1),1);
ind = max(globalnr_1h(:,4));

el5v = [ el4v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,5) = reshape(el5v,N*(N+1),1);
ind = max(globalnr_1v(:,5));

el5h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,5) = reshape(el5h',N*(N+1),1);
ind = max(globalnr_1h(:,5));

el6v = [ el5v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,6) = reshape(el6v,N*(N+1),1);
ind = max(globalnr_1v(:,6));

el6h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,6) = reshape(el6h',N*(N+1),1);
ind = max(globalnr_1h(:,6));

el7v = [ el6v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,7) = reshape(el7v,N*(N+1),1);
ind = max(globalnr_1v(:,7));

el7h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,7) = reshape(el7h',N*(N+1),1);
ind = max(globalnr_1h(:,7));

el8v = [ el7v(N+1,:) ; reshape(ind+(1:N*(N-1)),N-1,N) ; el1v(1,:) ];
globalnr_1v(:,8) = reshape(el8v,N*(N+1),1);
ind = max(globalnr_1v(:,8));

el8h = reshape(ind+(1:N*(N+1)),N+1,N)';
globalnr_1h(:,8) = reshape(el8h',N*(N+1),1);
ind = max(globalnr_1h(:,8));

el9v = reshape(ind+(1:N*(N+1)),N+1,N);
globalnr_1v(:,9) = reshape(el9v,N*(N+1),1);
ind = max(globalnr_1v(:,9));

el9h = [ el1h(:,N+1) reshape(ind+(1:N2),N,N)'  ];
globalnr_1h(:,9) = reshape(el9h',N*(N+1),1);
ind = max(globalnr_1h(:,9));

el10v = [ el9v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,10) = reshape(el10v,N*(N+1),1);
ind = max(globalnr_1v(:,10));

el10h = [ el2h(:,N+1) reshape(ind+(1:N2),N,N)'  ];
globalnr_1h(:,10) = reshape(el10h',N*(N+1),1);
ind = max(globalnr_1h(:,10));

el11v = [ el10v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,11) = reshape(el11v,N*(N+1),1);
ind = max(globalnr_1v(:,11));

el11h = [ el3h(:,N+1) reshape(ind+(1:N2),N,N)'  ];
globalnr_1h(:,11) = reshape(el11h',N*(N+1),1);
ind = max(globalnr_1h(:,11));

el12v = [ el11v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,12) = reshape(el12v,N*(N+1),1);
ind = max(globalnr_1v(:,12));

el12h = [ el4h(:,N+1) reshape(ind+(1:N2),N,N)'  ];
globalnr_1h(:,12) = reshape(el12h',N*(N+1),1);
ind = max(globalnr_1h(:,12));

el13v = [ el12v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,13) = reshape(el13v,N*(N+1),1);
ind = max(globalnr_1v(:,13));

el13h = [ el5h(:,N+1) reshape(ind+(1:N2),N,N)'  ];
globalnr_1h(:,13) = reshape(el13h',N*(N+1),1);
ind = max(globalnr_1h(:,13));

el14v = [ el13v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,14) = reshape(el14v,N*(N+1),1);
ind = max(globalnr_1v(:,14));

el14h = [ el6h(:,N+1) reshape(ind+(1:N2),N,N)'  ];
globalnr_1h(:,14) = reshape(el14h',N*(N+1),1);
ind = max(globalnr_1h(:,14));

el15v = [ el14v(N+1,:) ; reshape(ind+(1:N2),N,N) ];
globalnr_1v(:,15) = reshape(el15v,N*(N+1),1);
ind = max(globalnr_1v(:,15));

el15h = [ el7h(:,N+1) reshape(ind+(1:N2),N,N)'  ];
globalnr_1h(:,15) = reshape(el15h',N*(N+1),1);
ind = max(globalnr_1h(:,15));

el16v = [ el15v(N+1,:) ; reshape(ind+(1:N*(N-1)),N-1,N) ; el9v(1,:) ];
globalnr_1v(:,16) = reshape(el16v,N*(N+1),1);
ind = max(globalnr_1v(:,16));

el16h = [ el8h(:,N+1) reshape(ind+(1:N2),N,N)'  ];
globalnr_1h(:,16) = reshape(el16h',N*(N+1),1);
ind = max(globalnr_1h(:,16));


% 2-cells
globalnr_2 = zeros(N^2,16);
for i=1:16
    globalnr_2(:,i) = (i-1)*N2+(1:N2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_0 = max(max(globalnr_0));
nr_1 = max(max([globalnr_1v globalnr_1h]));
nr_2 = max(max(globalnr_2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%