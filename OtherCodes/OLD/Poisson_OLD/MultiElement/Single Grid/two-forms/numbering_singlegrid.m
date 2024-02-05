function [globalnr_0,globalnr_1h,globalnr_1v,globalnr_2] = numbering_singlegrid()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% [globalnr_0,globalnr_1h,globalnr_1v,globalnr_2]                         %
%                                      = numbering(N,numColumns,numRows)  %
%                                                                         %
% Global numbering for single grid problems in 2D                         %
%                                                                         %
% Global numbering for point values (0-cells)                             %
% Global numbering for lines values (1-cells)                             %
% Global numbering for surface values (2-cells)                           %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 28-10-2010                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N numRows numColumns

N2=N*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


globalnr_0 = zeros(numColumns*N+1,numRows*N+1);
% r==1
globalnr_0(1:N+1,1:N+1) = reshape(1:(N+1)^2,N+1,N+1);
localnr = reshape(1:N*(N+1),N,N+1);
for c=2:numColumns
    globalnr_0((c-1)*N+1+(1:N),1:N+1) = (c-1)*N*(N+1)+(N+1)+localnr;
end
localnr = reshape(1:N2,N,N);
for r=2:numRows
    for c=1:numColumns
        if c==1
            globalnr_0(1:N+1,(r-1)*N+1+(1:N)) = (N*numColumns+1)*(N*(r-1)+1) + reshape(1:N*(N+1),N+1,N);
        else
            globalnr_0((c-1)*N+1+(1:N),(r-1)*N+1+(1:N)) = (N*numColumns+1)*(N*(r-1)+1)+(c-1)*N2+N + localnr;
        end
    end
end

% disp('  '); disp('globalnr_0 = '); disp('  ');
% disp(num2str(flipud(globalnr_0')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalnr_1h = zeros(numColumns*N,numRows*N+1);
globalnr_1v = zeros(numColumns*N+1,numRows*N);
ind = 0;
for r=1:numRows
    for c=1:numColumns
        globalnr_1v((c-1)*N+(c>1)+(1:N+(c==1)),(r-1)*N+(1:N)) = ind + reshape(1:N*(N+(c==1)),N+(c==1),N) ;
        ind = ind + N2 + (c==1)*N;
        globalnr_1h((c-1)*N+(1:N),(r-1)*N+(r>1)+(1:N+(r==1))) = ind + reshape(1:N*(N+(r==1)),N+(r==1),N)';
        ind = ind + N*(N+(r==1));
    end
end

% disp('  '); disp('globalnr_1v = '); disp('  ');
% disp(num2str(flipud(globalnr_1v')));
% disp('  '); disp('globalnr_1h = '); disp('  ');
% disp(num2str(flipud(globalnr_1h')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localnr = reshape(1:N2,N,N);

globalnr_2 = zeros(numColumns*N,numRows*N);
for j=1:numRows
    for i=1:numColumns
        globalnr_2((i-1)*N+(1:N),(j-1)*N+(1:N)) = ((i-1)+(j-1)*numColumns)*N2+localnr;
    end
end

% disp('  '); disp('globalnr_2 = '); disp('  ');
% disp(num2str(flipud(globalnr_2')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%