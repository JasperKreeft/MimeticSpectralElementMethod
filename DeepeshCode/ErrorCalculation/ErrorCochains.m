function error = ErrorCochains(computedCochain, exactCochain)

% Error in cochains.
% error = ErrorCochains(computedCochain, exactCochain)
% Error = Sqrt(Sum(error.^2))
% where 
% error = (computedCochain - exactCochain)
%
% Inputs:
%
%
% computedCochain, exactCochain
%
%
% Outputs:
%
%     error
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $2/03/2012$

    % Reshape matrices to make sure they are the same size
    exactCochain = reshape(exactCochain, size(computedCochain));
    
    error = sum(sum((computedCochain - exactCochain).^2));

end