function []=numbering_CoVect(testcase)

% clear all
% clc
% testcase = 'square';
% N=1;
% numRows = 1;
% numColumns = 2;
% numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N numElements numRows numColumns
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global globalnr_11_xx globalnr_11_yx globalnr_11_yy globalnr_11_xy
global globalnr_21_x globalnr_21_y
global nr_0 nr_1 nr_2 nr_11 nr_21

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


% numbering 0,1-cells

% NOT CLEAR YET !!!!!



% numbering 1,1-cells
globalnr_11_xx = zeros(N*(N+2),numElements);
globalnr_11_yx = zeros((N+1)^2,numElements);
globalnr_11_yy = zeros(N*(N+2),numElements);
globalnr_11_xy = zeros((N+1)^2,numElements);

ind = 0;
inda = 1:N+2:N*(N+2);
indb = N+2:N+2:N*(N+2);
indc = 1:N*(N+2); indc(inda) = [];
indd = 1:N+1;
inde = N*(N+1)+1:(N+1)^2;
indf = 1:(N+1)^2; indf(indd) = [];
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        if c==1
globalnr_11_xx(:,rc) = ind + (1:N*(N+2));
ind = ind + N*(N+2);
        else
globalnr_11_xx(inda,rc) = globalnr_11_xx(indb,rc-1);
globalnr_11_xx(indc,rc) = ind+(1:N*(N+1));
ind = ind + N*(N+1);
        end
        
        if r==1
globalnr_11_yx(:,rc) = ind+(1:(N+1)^2);
ind = ind + (N+1)^2;
        else
globalnr_11_yx(indd,rc) = globalnr_11_yx(inde,rc-numColumns);
globalnr_11_yx(indf,rc) = ind+(1:N*(N+1));
ind = ind + N*(N+1);
        end
        
        if r==1
globalnr_11_yy(:,rc) = ind+(1:N*(N+2));
ind = ind + N*(N+2);
        else
globalnr_11_yy(inda,rc) = globalnr_11_yy(indb,rc-numColumns);
globalnr_11_yy(indc,rc) = ind+(1:N*(N+1));
ind = ind + N*(N+1);
        end
        
        if c==1
globalnr_11_xy(:,rc) = ind + (1:(N+1)^2);
ind = ind + (N+1)^2;
        else
globalnr_11_xy(indd,rc) = globalnr_11_xy(inde,rc-1);
globalnr_11_xy(indf,rc) = ind+(1:N*(N+1));
ind = ind + N*(N+1);
        end
        
    end
end


% numbering 2,1-cells
globalnr_21_x = zeros(N*(N+1),numElements);
globalnr_21_y = zeros(N*(N+1),numElements);
for i=1:numElements
    globalnr_21_x(:,i) = 2*(i-1)*N*(N+1)+(1:N*(N+1));
    globalnr_21_y(:,i) = (2*i-1)*N*(N+1)+(1:N*(N+1));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nr_0 = max(max(globalnr_0));
nr_1 = max(max([globalnr_1v globalnr_1h]));
nr_2 = max(max(globalnr_2));
nr_11 = max([ max(max(globalnr_11_xx)) max(max(globalnr_11_yx)) max(max(globalnr_11_yy)) max(max(globalnr_11_xy)) ]);
nr_21 = max(max([globalnr_21_x globalnr_21_y]));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%