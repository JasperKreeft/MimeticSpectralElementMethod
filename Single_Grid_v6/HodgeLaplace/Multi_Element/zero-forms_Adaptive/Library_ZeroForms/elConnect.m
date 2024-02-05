function [cnectN cnectL] = elConnect(hx,hy)
%
%   elConnect creates the element connectivity matrices in 2D
%
%   [cnectN cnectL] = elConnect(hx,hy)
%
%   input:
%       hx  :: number of elements in x-direction
%       hy  :: number of elements in y-direction
%
%   output:
%       cnectN  :: connectivity matrix with element numbers
%       cnectL  :: logical of cnectN
%
%   Copyright 2012 Peter Kuystermans
%   $Revision: 1.0 $  $Date: 18/03/2012 $  
%   Copyright 2012 Jasper Kreeft
%   $Revision: 1.1 $  $Date: 08/05/2012 $

%-------------------------------------------------------------------------%
% cnectN                                                                  %
%-------------------------------------------------------------------------%

    h = hx*hy;           % total number of elements
    cnectN = zeros(h,4); % creat empty cnectN

    % inner elements
    for i=2:hy-1
        for j=2:hx-1

            % element number 
            elNo = hx*(i-1)+j;        

            % store neighbours
            cnectN(elNo,:) = [ elNo-1 elNo+1 elNo-hx elNo+hx ]; % [ left right below above ]

        end
    end

    %%% domain corner boundary elements
    % lower left corner
    elNo = 1;
    if hx>1 % right connection (2)
        cnectN(elNo,2) = elNo+1;
    end
    if hy>1 % top connection (4)
        cnectN(elNo,4) = elNo+hx;
    end
    % lower right corner
    elNo = hx;
    if hy>1 % top connection (4)
        cnectN(elNo,4) = elNo+hx;
    end    
    if hx>1 % left connection (1)
        cnectN(elNo,1) = elNo-1;
    end
    % top right corner
    elNo = hx*hy;
    if hy>1 % bottom connection (3)
        cnectN(elNo,3) = elNo-hx;
    end 
    if hx>1 % left connection (1)
        cnectN(elNo,1) = elNo-1;
    end
    % top left corner
    elNo = hx*(hy-1)+1;
    if hy>1 % bottom connection (3)
        cnectN(elNo,3) = elNo-hx;
    end
    if hx>1 % right connection (2)
        cnectN(elNo,2) = elNo+1;
    end

    %%% domain non-corner/inner boundary elements
    % lower side
    elNo = (2:hx-1)';
    cnectN(elNo,[1 2]) = [elNo-1 elNo+1];
    if hy>1 % top connection (4)
        cnectN(elNo,4) = elNo+hx;
    end
    % right side
    elNo = (2*hx:hx:((hy-1)*hx))';
    cnectN(elNo,[3 4]) = [elNo-hx elNo+hx];
    if hx>1 % left connection (1)
        cnectN(elNo,1) = elNo-1;
    end
    % upper side
    elNo = ((hx*(hy-1)+2):(hx*hy-1))';
    cnectN(elNo,[1 2]) = [elNo-1 elNo+1];
    if hy>1 % bottom connection (3)
        cnectN(elNo,3) = elNo-hx;
    end
    % left side
    elNo = ((1+hx):hx:(hx*(hy-1)))';
    cnectN(elNo,[3 4]) = [elNo-hx elNo+hx];
    if hx>1 % right connection (2)
        cnectN(elNo,2) = elNo+1;
    end

%-------------------------------------------------------------------------%
% cnectLogical                                                            %
%-------------------------------------------------------------------------%

    cnectL = logical(cnectN>0);
    
end