% Marc's line polynomial
% histopolation

function [e]=EdgeVal(dhdx)

[ni nj] = size(dhdx);

e = zeros(ni-1,nj);
for i=1:ni-1
    for j=1:nj
        e(i,j) = -sum(dhdx(1:i,j));
    end
end