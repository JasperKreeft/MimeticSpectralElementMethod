function xMapped = MapParaToPhy(xParametric, xBound)

% xMapped = MapParaToPhy(xParametric, xBound)
%
% xParametric      :: points in parametric domain [nDomains X integrationPoints]
% xBound           :: domain boundaries [nDomains X 2]
%
% Takes as input the boundaries of a parametric domain (xBound), and on 
% each such domain, parametric points (xParametric) are mapped to a physical space
% (xMapped).
% Boundaries of parametric domain represent boundaries in xBound in physical
% space. xParametric can be anything between -1 and 1. Also, a
% linear map
% map := xPhy = c * xParametric + d
% is assumed.
%
% Copyright 2012 Deepesh Toshniwal
% Revision 1.0 $2/3/2012$

    %Make sure that matrices received are of correct dimensions
    xBound = reshape(xBound,size(xBound,1)*size(xBound,2)/2,2);
    xParametric = reshape(xParametric,size(xBound,1),size(xParametric,1)*size(xParametric,2)/size(xBound,1));

    %Linear Map (x = c(xi)+d)
    c = 0.5*(xBound(:,2)-xBound(:,1));
    d = 0.5*(xBound(:,2)+xBound(:,1));

    c = repmat(c,1,size(xParametric,2));
    d = repmat(d,1,size(xParametric,2));

    % generate points in physical space
    xMapped = c.*xParametric + d;
    
end