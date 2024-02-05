function uu = reconstruct_twoforms(U,e,J)

% only for single grid
% reconstruct 2-forms

global N numColumns numRows

if isempty(numColumns)
    numColumns = 1;
end
if isempty(numRows)
    numRows = 1;
end

nn = size(e,2);

uu = zeros(nn,nn,numColumns*numRows);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        uu(:,:,rc) = (e'*U((c-1)*N+(1:N),(r-1)*N+(1:N))*e).*J(:,:,rc);
    end
end
