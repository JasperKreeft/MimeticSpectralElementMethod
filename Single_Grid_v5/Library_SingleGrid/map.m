function [y,dy]=map(x,L,H)

q = (x*2/L);

% c = 15; % must be odd
% if even(c); disp('mapping coefficient changed'); c=c+1; end
% slope = 7/10;
% 
% y  = L/2*((1-slope)*q.^c+q*slope);
% dy = ((1-slope)*c*q.^(c-1)+slope)/H;

y  = q;
dy = ones(size(q))/H;