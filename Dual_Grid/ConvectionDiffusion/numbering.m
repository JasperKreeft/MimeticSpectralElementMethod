function [globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering()
% [globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering(N,numColumns,numRows)
% 
% Global numbering for point values (0-cells) on dual-grid
% Global numbering for lines values (1-cells) on primal-grid
% Global numbering for surface values (2-cells) on primal-grid

global N N2 numRows numColumns

N2=N*N;


localnr = reshape(1:N2,N,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_0 = zeros(numColumns*(N+1)+1,numRows*(N+1)+1);
for r=1:numRows
    for c=1:numColumns
        
        max_previous_element = (r-1)*(numColumns*(N2+2*N)+N)+...     % all previous rows
                               (r==1)*(c-1)*N+(r>1)*numColumns*N+... % bottom nodes
                               (c-1)*(N2+2*N)+...                    % elements at left-side                               
                               N*(c>1);                              % left boundary in row

        globalnr_0((c-1)*(N+1)+1+(1:N),(r-1)*(N+1)+1+(1:N)) = localnr+max_previous_element;

        if c==1
            globalnr_0(1,(r-1)*(N+1)+1+(1:N)) = (1:2:2*N)+N2+...
                          (r-1)*(numColumns*(N2+2*N)+N)+(r>1)*N*numColumns;

            globalnr_0(N+2,(r-1)*(N+1)+1+(1:N)) = (2:2:2*N)+N2+...
                          (r-1)*(numColumns*(N2+2*N)+N)+(r>1)*N*numColumns;
        else
            globalnr_0(c*(N+1)+1,(r-1)*(N+1)+1+(1:N)) = (1:N)+N2+(c-1)*(N2+2*N)+...
                          (r==1)*(c-1)*N+...
                          N*(c>1)+...
                          (r-1)*(numColumns*(N2+2*N)+N)+(r>1)*N*numColumns;
        end
        if r==1
            globalnr_0((c-1)*(N+1)+1+(1:N),1)   = (c-1)*(N2+3*N)+(N2+2*N)+(1:2:2*N);
            globalnr_0((c-1)*(N+1)+1+(1:N),N+2) = (c-1)*(N2+3*N)+(N2+2*N)+(2:2:2*N);
        else
            globalnr_0((c-1)*(N+1)+1+(1:N),r*(N+1)+1) = max_previous_element+N2+(1+(c==1))*N+(1:N);
        end
    end
end

% disp('  '); disp('globalnr_0 = '); disp('  ');
% disp(num2str(flipud(globalnr_0')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_2 = zeros(numColumns*N,numRows*N);
for j=1:numRows
    for i=1:numColumns
        globalnr_2((i-1)*N+(1:N),(j-1)*N+(1:N)) = ((i-1)+(j-1)*numColumns)*N2+localnr;
    end
end

% disp('  '); disp('globalnr_2 = '); disp('  '); disp(num2str(flipud(globalnr_2')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_1h = zeros(numColumns*(N+1),numRows*N);

localnr = reshape(1:N*(N+1),N+1,N);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        ind1 = (c-1)*(N+1)+(1:N+1);
        ind2 = (r-1)*N+(1:N);
        globalnr_1h(ind1,ind2) = localnr+(rc-1)*2*N*(N+1);
    end
end

% disp('  '); disp('globalnr_1h = '); disp('  ');
% disp(num2str(flipud(globalnr_1h')));



globalnr_1v = zeros(numColumns*N,numRows*(N+1));

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        ind1 = (c-1)*N+(1:N);
        ind2 = (r-1)*(N+1)+(1:N+1);
        globalnr_1v(ind1,ind2) = localnr'+(rc-1/2)*2*N*(N+1);
    end
end

% disp('  '); disp('globalnr_1v = '); disp('  ');
% disp(num2str(flipud(globalnr_1v')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%